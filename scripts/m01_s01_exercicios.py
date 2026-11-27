"""Modulo 1 / Sessao 1 - exercicios 2 a 6 da Aula 1, resolvidos.

Ex. 2  barcodes rbcL e matK -> FASTA (insumo do Modulo 3)
Ex. 3  quem depositou o registro (REFERENCE e comentario estruturado)
Ex. 4  copias anotadas: identicas, diferentes so na borda ou divergentes?
Ex. 5  comparacao congenerica: o que e biologia e o que e pratica de anotacao
Ex. 6  benchmark honesto: gc_janela O(n) contra a versao ingenua O(n * w / p)

Pre-requisito, uma vez: os registros no cache local.
    python scripts/m01_s01_plastoma.py
    python scripts/m01_s01_plastoma.py --accession NC_026450.1

Uso, a partir de 01_estudos com o .venv ativo:
    python scripts/m01_s01_exercicios.py
    python scripts/m01_s01_exercicios.py --sem-benchmark
    python scripts/m01_s01_exercicios.py --focal KY404932.1 --comparar NC_026450.1

Saidas em BIOINFO_DADOS (03_dados)/processados/plastomas:
    <focal>_barcodes.fasta
    <accession>_copias_pareadas.parquet   (um por registro)
    m01_s01_benchmark_gc_janela.parquet
"""

from __future__ import annotations

# bioflora antes de qualquer pacote cientifico (teto de threads do BLAS).
# O comentario "isort: split" impede o ruff de mover este import para depois
# do pandas, o que reativaria o [AVISO] de threads.
import bioflora  # noqa: F401

# isort: split
import argparse
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from bioflora.config import dir_dados
from bioflora.ncbi import caminho_cache
from bioflora.recursos import cronometrar, log, log_memoria
from bioflora.sequencias import (
    comparar_copias,
    gc_janela,
    gc_janela_ingenua,
    genes_repetidos,
    qualificador,
    resumo_registro,
)

BARCODES = ("rbcL", "matK")
JANELA = 1000
PASSOS = (500, 250, 100, 50, 10, 1)  # w/p de 2 a 1.000


def titulo(texto: str) -> None:
    print()
    print("=" * 72)
    print(texto)
    print("=" * 72)


def ler_registro(accession: str) -> SeqRecord | None:
    """GenBank do cache local ou None, com a instrucao para baixar."""
    arq = caminho_cache(accession)
    if not arq.exists():
        log("ERRO", f"{arq} nao existe")
        log("ERRO", f"rode: python scripts/m01_s01_plastoma.py --accession {accession}")
        return None
    return SeqIO.read(arq, "genbank")


def nomes_de_genes(record: SeqRecord, tipo: str) -> set[str]:
    return {
        qualificador(f, "gene")
        for f in record.features
        if f.type == tipo and qualificador(f, "gene")
    }


# --- Exercicio 2 -------------------------------------------------------------------


def ex2_barcodes(record: SeqRecord, saida: Path) -> None:
    titulo("EXERCICIO 2: barcodes rbcL e matK -> FASTA")
    organismo = record.annotations.get("organism", record.id)
    canonicos = {nome.lower(): nome for nome in BARCODES}
    registros = []
    for f in record.features:
        anotado = qualificador(f, "gene")
        gene = canonicos.get(anotado.lower())
        if f.type != "CDS" or gene is None:
            continue
        if anotado != gene:  # ex.: "rbcl" em NC_026450.1
            log("AVISO", f"{gene} anotado como '{anotado}' (nome fora do padrao)")
        nt = f.extract(record.seq)
        n = len(nt)
        nivel = "OK" if n % 3 == 0 else "AVISO"
        extremos = f"inicio {nt[:3]} | fim {nt[-3:]}"
        log(nivel, f"{gene}: {n} pb | {n // 3} codons | {extremos}")
        registros.append(
            SeqRecord(
                nt,
                id=f"{record.id}_{gene}",
                description=f"{organismo} {gene} CDS ({n} pb)",
            )
        )
    encontrados = {r.id.rsplit("_", 1)[1] for r in registros}
    faltam = sorted(set(BARCODES) - encontrados)
    if faltam:
        log("AVISO", f"sem CDS anotada para: {', '.join(faltam)}")
    arq = saida / f"{record.id}_barcodes.fasta"
    SeqIO.write(registros, arq, "fasta")
    log("OK", f"{len(registros)} sequencia(s) -> {arq}")


# --- Exercicio 3 -------------------------------------------------------------------


def ex3_deposito(record: SeqRecord) -> None:
    titulo(f"EXERCICIO 3: quem depositou {record.id}?")
    for i, ref in enumerate(record.annotations.get("references", []), start=1):
        autores = ref.authors or ref.consrtm or "-"
        log("OK", f"referencia {i}: {autores}")
        print(f"     titulo : {ref.title}")
        print(f"     revista: {ref.journal}")
        if ref.pubmed_id:
            print(f"     PubMed : {ref.pubmed_id}")

    # Registros mais novos trazem o metodo em ##Assembly-Data-START## no COMMENT.
    estruturado = record.annotations.get("structured_comment", {})
    montagem = estruturado.get("Assembly-Data", {})
    if montagem:
        for chave, valor in montagem.items():
            log("OK", f"{chave}: {valor}")
    else:
        log("AVISO", "sem comentario estruturado de montagem (Assembly-Data):")
        log("AVISO", "o registro nao declara tecnologia nem metodo de sequenciamento")
    comentario = " ".join(record.annotations.get("comment", "").split())
    if comentario:
        log("OK", f"COMMENT: {comentario[:160]}")


# --- Exercicio 4 -------------------------------------------------------------------


def ex4_copias(record: SeqRecord, saida: Path) -> pd.DataFrame:
    titulo(f"EXERCICIO 4: copias anotadas de {record.id}")
    pares = comparar_copias(record)
    if pares.empty:
        log("OK", "nenhum gene anotado em mais de uma copia")
        return pares
    for p in pares.itertuples(index=False):
        log(
            "OK",
            f"{p.gene:<9} {p.copia_a}x{p.copia_b} | {p.orientacao:<9} | "
            f"{p.relacao:<10} | {p.comprimento_a} x {p.comprimento_b} pb | "
            f"{p.substituicoes} subst. | {p.gaps} gaps | {p.identidade_pct:.1f}%",
        )
    log("OK", f"pares por relacao: {pares['relacao'].value_counts().to_dict()}")
    arq = saida / f"{record.id}_copias_pareadas.parquet"
    pares.to_parquet(arq, index=False)
    log("OK", f"tabela de pares -> {arq}")
    return pares


# --- Exercicio 5 -------------------------------------------------------------------


def ex5_comparacao(records: list[SeqRecord]) -> pd.DataFrame:
    titulo("EXERCICIO 5: comparacao congenerica (biologia ou anotacao?)")
    linhas = []
    for rec in records:
        r = resumo_registro(rec)
        rep = genes_repetidos(rec)
        linhas.append(
            {
                "accession": rec.id,
                "especie": r["organismo"],
                "pb": r["comprimento_pb"],
                "gc_pct": r["gc_pct"],
                "cds": r["genes_unicos"]["CDS"],
                "trna": r["genes_unicos"]["tRNA"],
                "rrna": r["genes_unicos"]["rRNA"],
                "feat_exon": r["n_features"].get("exon", 0),
                "feat_intron": r["n_features"].get("intron", 0),
                "repetidos": ", ".join(
                    f"{g} {n}x" for _, g, n in rep.itertuples(index=False)
                ),
            }
        )
    tabela = pd.DataFrame(linhas)
    print(tabela.drop(columns="repetidos").to_string(index=False))
    for linha in linhas:
        log("OK", f"{linha['accession']} repetidos: {linha['repetidos'] or '-'}")

    a, b = records[0], records[1]
    for tipo in ("CDS", "tRNA", "rRNA"):
        na, nb = nomes_de_genes(a, tipo), nomes_de_genes(b, tipo)
        so_a, so_b = na - nb, nb - na
        # Nomes que so diferem na caixa (rbcL x rbcl) sao o mesmo gene.
        caixa = {x for x in so_a if x.lower() in {y.lower() for y in so_b}}
        so_a -= caixa
        so_b = {y for y in so_b if y.lower() not in {x.lower() for x in caixa}}
        if caixa:
            log("AVISO", f"{tipo}: nomes que so diferem na caixa: {sorted(caixa)}")
        if so_a or so_b:
            log(
                "AVISO",
                f"{tipo}: so em {a.id}: {', '.join(sorted(so_a)) or '-'} | "
                f"so em {b.id}: {', '.join(sorted(so_b)) or '-'}",
            )
        else:
            log("OK", f"{tipo}: mesmo conjunto de genes nos dois registros")
    log("OK", "features 'exon' x 'intron' diferentes = estilos de anotacao diferentes")
    log("OK", "nome diferente nao prova gene diferente: a Sessao 2 compara sequencias")
    return tabela


# --- Exercicio 6 -------------------------------------------------------------------


def ex6_benchmark(record: SeqRecord, saida: Path) -> pd.DataFrame | None:
    titulo("EXERCICIO 6: gc_janela O(n) contra a ingenua O(n * w / p)")
    seq = record.seq
    log("OK", f"{len(seq):,} pb | janela {JANELA} | passos {PASSOS}")
    log("OK", "cronometrando (pode levar ate 1 min no notebook)...")
    linhas = []
    for passo in PASSOS:
        # 1. Antes de cronometrar, provar que as duas respondem a mesma coisa.
        rapida = gc_janela(seq, JANELA, passo)
        ingenua = gc_janela_ingenua(seq, JANELA, passo)
        mesmos_centros = np.array_equal(rapida[0], ingenua[0])
        mesmo_gc = np.allclose(rapida[1], ingenua[1], equal_nan=True)
        if not (mesmos_centros and mesmo_gc):
            log("ERRO", f"passo {passo}: resultados diferentes; benchmark invalido")
            return None
        # 2. Cronometrar: menor de 3 lotes, lote dimensionado por autorange.
        ms_prefixo = cronometrar(lambda p=passo: gc_janela(seq, JANELA, p), 3)
        ms_ingenua = cronometrar(lambda p=passo: gc_janela_ingenua(seq, JANELA, p), 3)
        linhas.append(
            {
                "janela": JANELA,
                "passo": passo,
                "w_sobre_p": JANELA / passo,
                "n_janelas": int(rapida[0].size),
                "ms_prefixo": ms_prefixo,
                "ms_ingenua": ms_ingenua,
                "razao": ms_ingenua / ms_prefixo,
            }
        )
    tab = pd.DataFrame(linhas)
    tab["accession"] = record.id
    tab["python"] = platform.python_version()
    tab["numpy"] = np.__version__
    tab["maquina"] = f"{platform.system()} {platform.machine()}"
    print(
        tab[["passo", "w_sobre_p", "n_janelas", "ms_prefixo", "ms_ingenua", "razao"]]
        .round(3)
        .to_string(index=False)
    )

    vence = tab[tab["razao"] > 1]
    if vence.empty:
        log("AVISO", "a versao O(n) nao venceu em nenhum passo testado")
    else:
        primeiro = vence.iloc[0]
        log(
            "OK",
            f"a versao O(n) vence a partir de w/p = {primeiro['w_sobre_p']:.0f} "
            f"(passo {int(primeiro['passo'])}); em w/p = 1.000 e "
            f"{tab['razao'].iloc[-1]:.0f}x mais rapida",
        )
    variacao = tab["ms_prefixo"].max() / tab["ms_prefixo"].min()
    log("OK", f"custo da O(n) varia {variacao:.1f}x quando w/p varia 500x")
    arq = saida / "m01_s01_benchmark_gc_janela.parquet"
    tab.to_parquet(arq, index=False)
    log("OK", f"benchmark -> {arq}")
    return tab


def main() -> int:
    parser = argparse.ArgumentParser(description="M1/S1 - exercicios 2 a 6 da Aula 1")
    parser.add_argument("--focal", default="KY404932.1")
    parser.add_argument("--comparar", nargs="*", default=["NC_026450.1"])
    parser.add_argument("--sem-benchmark", action="store_true")
    args = parser.parse_args()

    focal = ler_registro(args.focal)
    if focal is None:
        return 1
    outros = [ler_registro(acc) for acc in args.comparar]
    faltando = any(r is None for r in outros)
    outros = [r for r in outros if r is not None]

    saida = dir_dados() / "processados" / "plastomas"
    saida.mkdir(parents=True, exist_ok=True)

    ex2_barcodes(focal, saida)
    ex3_deposito(focal)
    for rec in [focal, *outros]:
        ex4_copias(rec, saida)
    if outros:
        ex5_comparacao([focal, *outros])
    else:
        log("AVISO", "exercicio 5 pulado: nenhum registro para comparar")
    if not args.sem_benchmark:
        ex6_benchmark(focal, saida)

    log_memoria("RAM ao final")
    return 1 if faltando else 0


if __name__ == "__main__":
    sys.exit(main())
