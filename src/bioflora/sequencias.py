"""Sequencias como dados: coordenadas, composicao e anotacoes (Modulo 1).

Convencoes que este modulo deixa explicitas:

- GenBank usa base 1 e intervalo FECHADO: ``[a..b]``.
- Biopython (e Python) usa base 0 e intervalo SEMIABERTO: ``[a-1, b)``.
- ``len(feature.location)`` soma as partes (exons) de uma localizacao composta.
- ``feature.extract(seq)`` ja devolve a sequencia na orientacao 5'->3' do gene,
  inclusive para genes na fita negativa e genes com introns.
"""

from __future__ import annotations

from collections import Counter
from typing import Any

import numpy as np
import pandas as pd
from Bio.Data.CodonTable import TranslationError
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

TIPOS_PADRAO: tuple[str, ...] = ("gene", "CDS", "tRNA", "rRNA")
TABELA_PLASTIDEO = 11  # NCBI: "Bacterial, Archaeal and Plant Plastid Code"

_BASES_VALIDAS = np.frombuffer(b"ACGT", dtype=np.uint8)
_GC = np.frombuffer(b"GC", dtype=np.uint8)


# --- Composicao ------------------------------------------------------------------


def gc_percent(seq: Any) -> float:
    """Conteudo GC em % ignorando simbolos ambiguos (N, R, Y...)."""
    return 100.0 * gc_fraction(seq, ambiguous="remove")


def gc_janela(
    seq: Any, janela: int = 1000, passo: int = 500
) -> tuple[np.ndarray, np.ndarray]:
    """GC (%) em janelas deslizantes, em O(n), via somas de prefixo.

    Para a contagem acumulada C (C[0] = 0), o numero de G/C na janela [i, i+w) e
    C[i+w] - C[i]. O denominador conta so A/C/G/T, entao N nao dilui o GC.
    Devolve (centros das janelas em pb, base 1; GC em %).
    """
    if janela <= 0 or passo <= 0:
        raise ValueError("janela e passo devem ser positivos")
    dados = seq.encode("ascii") if isinstance(seq, str) else bytes(seq)
    arr = np.frombuffer(dados.upper(), dtype=np.uint8)
    if arr.size < janela:
        return np.array([], dtype=float), np.array([], dtype=float)

    eh_gc = np.isin(arr, _GC)
    valido = np.isin(arr, _BASES_VALIDAS)
    acum_gc = np.concatenate(([0], np.cumsum(eh_gc, dtype=np.int64)))
    acum_val = np.concatenate(([0], np.cumsum(valido, dtype=np.int64)))

    inicios = np.arange(0, arr.size - janela + 1, passo)
    n_gc = acum_gc[inicios + janela] - acum_gc[inicios]
    n_val = acum_val[inicios + janela] - acum_val[inicios]
    pct = np.full(inicios.size, np.nan)
    np.divide(100.0 * n_gc, n_val, out=pct, where=n_val > 0)
    centros = inicios + janela / 2 + 0.5
    return centros, pct


# --- Coordenadas e anotacoes ----------------------------------------------------


def coordenadas_genbank(feature: SeqFeature) -> tuple[int, int]:
    """(inicio, fim) como o GenBank mostra: base 1, intervalo fechado."""
    return int(feature.location.start) + 1, int(feature.location.end)


def qualificador(feature: SeqFeature, chave: str, padrao: str = "") -> str:
    """Primeiro valor de um qualificador da anotacao (ex.: /gene, /product)."""
    return feature.qualifiers.get(chave, [padrao])[0]


def features_para_dataframe(
    record: SeqRecord, tipos: tuple[str, ...] = TIPOS_PADRAO
) -> pd.DataFrame:
    """Tabela de anotacoes, uma linha por feature dos tipos pedidos.

    ``fita`` vale 1, -1 ou 0; 0 indica partes em fitas diferentes (ex.: genes com
    trans-splicing, como o rps12 de muitos plastomas).
    """
    linhas = []
    for f in record.features:
        if f.type not in tipos:
            continue
        inicio, fim = coordenadas_genbank(f)
        trecho = f.extract(record.seq)
        linhas.append(
            {
                "tipo": f.type,
                "gene": qualificador(f, "gene"),
                "produto": qualificador(f, "product"),
                "locus_tag": qualificador(f, "locus_tag"),
                "inicio": inicio,
                "fim": fim,
                "fita": int(f.location.strand or 0),
                "comprimento": len(f.location),
                "n_partes": len(f.location.parts),
                "gc_pct": round(gc_percent(trecho), 2) if len(trecho) else np.nan,
                "localizacao": str(f.location),
            }
        )
    colunas = [
        "tipo", "gene", "produto", "locus_tag", "inicio", "fim", "fita",
        "comprimento", "n_partes", "gc_pct", "localizacao",
    ]  # fmt: skip
    df = pd.DataFrame(linhas, columns=colunas)
    return df.astype(
        {
            "inicio": "int64",
            "fim": "int64",
            "fita": "int8",
            "comprimento": "int64",
            "n_partes": "int16",
            "gc_pct": "float64",
        }
    )


def genes_repetidos(
    record: SeqRecord, tipos: tuple[str, ...] = ("CDS", "tRNA", "rRNA")
) -> pd.DataFrame:
    """Genes anotados mais de uma vez por tipo (copias), como DataFrame ordenado."""
    contagem = Counter(
        (f.type, qualificador(f, "gene"))
        for f in record.features
        if f.type in tipos and qualificador(f, "gene")
    )
    linhas = [
        {"tipo": t, "gene": g, "copias": n} for (t, g), n in contagem.items() if n > 1
    ]
    df = pd.DataFrame(linhas, columns=["tipo", "gene", "copias"])
    return df.sort_values(["tipo", "gene"], ignore_index=True)


def resumo_registro(record: SeqRecord) -> dict[str, Any]:
    """Metadados e composicao basica de um registro GenBank."""
    contagem = Counter(f.type for f in record.features)
    genes_unicos = {
        tipo: len(
            {
                qualificador(f, "gene")
                for f in record.features
                if f.type == tipo and qualificador(f, "gene")
            }
        )
        for tipo in ("CDS", "tRNA", "rRNA")
    }
    return {
        "accession": record.id,
        "descricao": record.description,
        "organismo": record.annotations.get("organism", ""),
        "taxonomia": "; ".join(record.annotations.get("taxonomy", [])),
        "molecula": record.annotations.get("molecule_type", ""),
        "topologia": record.annotations.get("topology", ""),
        "comprimento_pb": len(record.seq),
        "gc_pct": round(gc_percent(record.seq), 2),
        "n_features": dict(contagem),
        "genes_unicos": genes_unicos,
    }


# --- Traducao e validacao de CDS ------------------------------------------------


def _diferencas(a: str, b: str) -> int | None:
    """Numero de posicoes diferentes (Hamming) se os comprimentos forem iguais."""
    if len(a) != len(b):
        return None
    return sum(x != y for x, y in zip(a, b, strict=True))


def validar_cds(record: SeqRecord, tabela: int = TABELA_PLASTIDEO) -> pd.DataFrame:
    """Traduz cada CDS a partir do DNA e confronta com a anotacao /translation.

    ``cds_valida`` usa ``translate(cds=True)``: exige codon de inicio valido para a
    tabela, comprimento multiplo de 3, um unico stop no fim e nenhum stop interno.
    Diferencas de poucos aminoacidos contra a anotacao sao a assinatura classica de
    edicao de RNA (C->U) em plastideos: o DNA nao e a proteina.
    Respeita /codon_start e /transl_table quando presentes.
    """
    linhas = []
    for f in record.features:
        if f.type != "CDS":
            continue
        tabela_f = int(qualificador(f, "transl_table", str(tabela)))
        inicio_leitura = int(qualificador(f, "codon_start", "1")) - 1
        nt = f.extract(record.seq)[inicio_leitura:]
        nt3 = nt[: len(nt) - len(nt) % 3]
        anotada = qualificador(f, "translation")

        try:
            proteina = str(nt.translate(table=tabela_f, cds=True))
            valida, motivo = True, ""
        except TranslationError as erro:
            proteina = str(nt3.translate(table=tabela_f)).rstrip("*")
            valida, motivo = False, str(erro)

        traduzida_bruta = str(nt3.translate(table=tabela_f))
        stops_internos = traduzida_bruta.rstrip("*").count("*")
        linhas.append(
            {
                "gene": qualificador(f, "gene"),
                "comprimento_nt": len(nt),
                "multiplo_de_3": len(nt) % 3 == 0,
                "codon_inicio": str(nt[:3]).upper(),
                "cds_valida": valida,
                "motivo": motivo,
                "stops_internos": stops_internos,
                "tem_anotacao": bool(anotada),
                "confere_anotacao": (proteina == anotada) if anotada else pd.NA,
                "diferencas_aa": _diferencas(proteina, anotada) if anotada else None,
            }
        )
    colunas = [
        "gene", "comprimento_nt", "multiplo_de_3", "codon_inicio", "cds_valida",
        "motivo", "stops_internos", "tem_anotacao", "confere_anotacao", "diferencas_aa",
    ]  # fmt: skip
    df = pd.DataFrame(linhas, columns=colunas)
    return df.astype(
        {
            "comprimento_nt": "int64",
            "stops_internos": "int64",
            "confere_anotacao": "boolean",
            "diferencas_aa": "Int64",
        }
    )
