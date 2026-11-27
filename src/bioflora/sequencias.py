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


def _tabela_de_consulta(simbolos: bytes) -> np.ndarray:
    """Tabela de 256 posicoes: 1 nos bytes de ``simbolos`` (maiusculos e
    minusculos), 0 no resto. Indexar por byte troca uma busca por um acesso direto.
    """
    tabela = np.zeros(256, dtype=np.uint8)
    codigos = np.frombuffer(simbolos.upper() + simbolos.lower(), dtype=np.uint8)
    tabela[codigos] = 1
    return tabela


_LUT_GC = _tabela_de_consulta(b"GC")
_LUT_VALIDAS = _tabela_de_consulta(b"ACGT")


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
    Maiusculas e minusculas valem igual. Devolve (centros das janelas em pb,
    base 1; GC em %).

    Custo: tres passagens lineares sobre a sequencia (consulta por byte e duas
    somas acumuladas), independentemente de ``janela`` e ``passo``. A versao
    ingenua custa O(n * janela / passo); por isso ela so perde quando
    janela / passo e grande (ver o exercicio 6 da Aula 1).
    """
    if janela <= 0 or passo <= 0:
        raise ValueError("janela e passo devem ser positivos")
    dados = seq.encode("ascii") if isinstance(seq, str) else bytes(seq)
    arr = np.frombuffer(dados, dtype=np.uint8)
    if arr.size < janela:
        return np.array([], dtype=float), np.array([], dtype=float)

    # int32 basta ate ~2,1 Gpb; acima disso, int64 evita estouro da soma.
    tipo_acum = np.int32 if arr.size < np.iinfo(np.int32).max else np.int64
    acum_gc = np.zeros(arr.size + 1, dtype=tipo_acum)
    acum_val = np.zeros(arr.size + 1, dtype=tipo_acum)
    np.cumsum(np.take(_LUT_GC, arr), dtype=tipo_acum, out=acum_gc[1:])
    np.cumsum(np.take(_LUT_VALIDAS, arr), dtype=tipo_acum, out=acum_val[1:])

    inicios = np.arange(0, arr.size - janela + 1, passo)
    n_gc = acum_gc[inicios + janela] - acum_gc[inicios]
    n_val = acum_val[inicios + janela] - acum_val[inicios]
    pct = np.full(inicios.size, np.nan)
    np.divide(100.0 * n_gc, n_val, out=pct, where=n_val > 0)
    centros = inicios + janela / 2 + 0.5
    return centros, pct


def gc_janela_ingenua(
    seq: Any, janela: int = 1000, passo: int = 500
) -> tuple[np.ndarray, np.ndarray]:
    """Referencia em O(n * janela / passo): reconta cada janela do zero.

    Mesma semantica de :func:`gc_janela` (centros em base 1, denominador A/C/G/T,
    caixa ignorada), para servir de oraculo nos testes e de linha de base em
    benchmarks. Cada janela custa ``janela`` leituras, mas em C (``str.count``),
    e as janelas sao (n - janela) / passo + 1. Nao use em producao.
    """
    if janela <= 0 or passo <= 0:
        raise ValueError("janela e passo devem ser positivos")
    texto = str(seq).upper()
    centros: list[float] = []
    pct: list[float] = []
    for inicio in range(0, len(texto) - janela + 1, passo):
        trecho = texto[inicio : inicio + janela]
        n_gc = trecho.count("G") + trecho.count("C")
        n_val = n_gc + trecho.count("A") + trecho.count("T")
        centros.append(inicio + janela / 2 + 0.5)
        pct.append(100.0 * n_gc / n_val if n_val else np.nan)
    return np.array(centros, dtype=float), np.array(pct, dtype=float)


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


def _alinhador_global():
    """Alinhamento global simples para comparar copias curtas (tRNA, rRNA, CDS)."""
    from Bio.Align import PairwiseAligner

    return PairwiseAligner(
        mode="global",
        match_score=1,
        mismatch_score=-1,
        open_gap_score=-2,
        extend_gap_score=-0.5,
    )


def comparar_copias(
    record: SeqRecord, tipos: tuple[str, ...] = ("CDS", "tRNA", "rRNA")
) -> pd.DataFrame:
    """Compara, par a par, as copias de cada gene anotado mais de uma vez.

    Responde a pergunta que a igualdade de strings nao responde: copias com
    "sequencias distintas" diferem na biologia ou so nos limites da anotacao?

    ``relacao``: ``identica`` (mesma sequencia), ``borda`` (uma copia contida na
    outra: diferem so nas pontas, tipico de limites anotados de forma diferente) ou
    ``divergente`` (substituicoes ou indels internos).
    ``orientacao``: ``direta`` (mesma fita) ou ``invertida`` (fitas opostas). As
    sequencias vem de ``feature.extract``, ja no sentido 5'->3' do gene, entao a
    comparacao vale para qualquer orientacao.
    ``identidade_pct`` = identidades / colunas do alinhamento global, incluindo
    gaps das pontas (por isso uma copia ``borda`` fica abaixo de 100%).
    """
    copias: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for f in record.features:
        gene = qualificador(f, "gene")
        if f.type not in tipos or not gene:
            continue
        inicio, fim = coordenadas_genbank(f)
        copias.setdefault((f.type, gene), []).append(
            {
                "inicio": inicio,
                "fim": fim,
                "fita": int(f.location.strand or 0),
                "seq": str(f.extract(record.seq)).upper(),
            }
        )

    alinhador = _alinhador_global()
    linhas = []
    for (tipo, gene), lista in sorted(copias.items()):
        lista = sorted(lista, key=lambda c: c["inicio"])
        for i in range(len(lista)):
            for j in range(i + 1, len(lista)):
                a, b = lista[i], lista[j]
                sa, sb = a["seq"], b["seq"]
                # So os pares divergentes precisam de alinhamento (custo O(la * lb));
                # nos outros dois casos a contagem sai direto dos comprimentos.
                if sa == sb:
                    relacao, ident, subst, gaps = "identica", len(sa), 0, 0
                elif sa in sb or sb in sa:
                    relacao = "borda"
                    ident, subst = min(len(sa), len(sb)), 0
                    gaps = abs(len(sa) - len(sb))
                else:
                    relacao = "divergente"
                    contagem = alinhador.align(sa, sb)[0].counts()
                    ident = contagem.identities
                    subst = contagem.mismatches
                    gaps = contagem.gaps
                colunas = ident + subst + gaps
                if 0 in (a["fita"], b["fita"]):
                    orientacao = "indefinida"
                else:
                    orientacao = "direta" if a["fita"] == b["fita"] else "invertida"
                linhas.append(
                    {
                        "tipo": tipo,
                        "gene": gene,
                        "copia_a": i + 1,
                        "copia_b": j + 1,
                        "inicio_a": a["inicio"],
                        "inicio_b": b["inicio"],
                        "fita_a": a["fita"],
                        "fita_b": b["fita"],
                        "comprimento_a": len(sa),
                        "comprimento_b": len(sb),
                        "orientacao": orientacao,
                        "relacao": relacao,
                        "identidades": ident,
                        "substituicoes": subst,
                        "gaps": gaps,
                        "identidade_pct": (
                            round(100 * ident / colunas, 2) if colunas else np.nan
                        ),
                    }
                )
    colunas_df = [
        "tipo", "gene", "copia_a", "copia_b", "inicio_a", "inicio_b", "fita_a",
        "fita_b", "comprimento_a", "comprimento_b", "orientacao", "relacao",
        "identidades", "substituicoes", "gaps", "identidade_pct",
    ]  # fmt: skip
    return pd.DataFrame(linhas, columns=colunas_df)


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
