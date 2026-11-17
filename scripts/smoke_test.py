"""Teste de fumaca do ambiente Python (Aula 0).

Uso, a partir de 01_estudos com o .venv ativo:
    python scripts/smoke_test.py            # inclui teste online (NCBI Entrez)
    python scripts/smoke_test.py --offline  # sem internet

Sai com codigo 1 se algum item critico falhar.
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys

# bioflora antes de qualquer pacote cientifico: aplica o teto de threads antes de o
# numpy ser carregado (pandas e Bio so sao importados dentro das funcoes abaixo)
import bioflora
from bioflora.config import RAIZ, configurar_entrez, dir_dados, dir_tmp
from bioflora.recursos import conectar_duckdb, log, log_memoria

ERROS: list[str] = []


def falha(mensagem: str) -> None:
    log("ERRO", mensagem)
    ERROS.append(mensagem)


def checar_sessao() -> None:
    log("OK", f"Python {sys.version.split()[0]} em {sys.executable}")
    if sys.prefix != sys.base_prefix:
        log("OK", f"venv ativo: {sys.prefix}")
    elif os.environ.get("CI"):
        # GitHub Actions instala direto no Python do runner, que ja e isolado
        log("AVISO", "fora de um venv (aceito em CI)")
    else:
        falha("fora de um venv: ative o .venv ou use .venv\\Scripts\\python.exe")
    threads = os.environ.get("OMP_NUM_THREADS")
    log("OK", f"bioflora {bioflora.__version__} | OMP_NUM_THREADS={threads}")

    # No workspace, o .env aponta para 03_dados; o fallback ./dados e so para clones
    env_arq = RAIZ / ".env"
    if env_arq.exists():
        log("OK", f".env encontrado: {env_arq}")
    else:
        log("AVISO", f".env ausente em {env_arq}: rode o Passo 2 de novo")
    for var, destino in (("BIOINFO_DADOS", dir_dados()), ("BIOINFO_TMP", dir_tmp())):
        if os.environ.get(var):
            log("OK", f"{var:<13} -> {destino}")
        else:
            log("AVISO", f"{var:<13} -> {destino} (fallback; falta no .env)")


def checar_pacotes() -> None:
    for nome in ("Bio", "numpy", "pandas", "pyarrow", "duckdb", "matplotlib", "psutil"):
        try:
            mod = importlib.import_module(nome)
            log("OK", f"{nome:<11} {getattr(mod, '__version__', '?')}")
        except Exception as exc:
            falha(f"{nome} nao importa: {exc}")


def checar_biopython() -> None:
    from Bio.Seq import Seq
    from Bio.SeqUtils import gc_fraction

    # Trecho no estilo do inicio do rbcL: codifica MSPQTETKASVGFKAGVKDY, o N-terminal
    # classico da subunidade grande da RuBisCO. Tabela 11 = codigo de plastideo.
    seq = Seq("ATGTCACCACAAACAGAGACTAAAGCAAGTGTTGGATTCAAAGCTGGTGTTAAAGATTAC")
    proteina = seq.translate(table=11)
    rc = seq.reverse_complement()
    gc = gc_fraction(seq)
    if str(proteina).startswith("MSPQTETKA") and len(rc) == len(seq):
        log("OK", f"Biopython: traducao tabela 11 = {proteina} | GC = {gc:.3f}")
    else:
        falha(f"Biopython: traducao inesperada {proteina}")


def checar_duckdb() -> None:
    import pandas as pd

    con = conectar_duckdb(memoria="2GB")
    try:
        destino = (dir_tmp() / "smoke.parquet").as_posix()
        df = pd.DataFrame(
            {"especie": ["Araucaria angustifolia", "Dalbergia nigra"], "n": [10, 20]}
        )
        con.register("df", df)
        con.execute(f"COPY df TO '{destino}' (FORMAT parquet)")
        total = con.execute(f"SELECT sum(n) FROM '{destino}'").fetchone()[0]
        limite = con.execute("SELECT current_setting('memory_limit')").fetchone()[0]
        threads = con.execute("SELECT current_setting('threads')").fetchone()[0]
        if total == 30:
            log("OK", f"DuckDB: Parquet ok | memory_limit={limite} | threads={threads}")
        else:
            falha(f"DuckDB: soma inesperada {total}")
    finally:
        con.close()


def checar_entrez(offline: bool) -> None:
    if offline:
        log("AVISO", "teste online do NCBI pulado (--offline)")
        return
    if not configurar_entrez():
        log("AVISO", "NCBI_EMAIL vazio no .env; Entrez pulado (preencha antes do M1)")
        return
    from Bio import Entrez

    termo = (
        '"Araucaria angustifolia"[Organism] AND chloroplast[Title] '
        'AND "complete genome"[Title]'
    )
    try:
        with Entrez.esearch(db="nucleotide", term=termo, retmax=3) as h:
            res = Entrez.read(h)
        ids = res["IdList"]
        log("OK", f"Entrez esearch: {res['Count']} plastoma(s) de A. angustifolia")
        if ids:
            with Entrez.esummary(db="nucleotide", id=ids[0]) as h:
                resumo = Entrez.read(h)[0]
            pb = int(resumo["Length"])
            log("OK", f"Entrez esummary: {resumo['AccessionVersion']} | {pb:,} pb")
    except Exception as exc:
        log("AVISO", f"NCBI indisponivel agora (nao critico): {exc}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Teste de fumaca do ambiente Python")
    parser.add_argument("--offline", action="store_true", help="sem internet")
    args = parser.parse_args()

    for etapa in (checar_sessao, checar_pacotes, checar_biopython, checar_duckdb):
        try:
            etapa()
        except Exception as exc:
            falha(f"{etapa.__name__}: {exc}")
    checar_entrez(args.offline)
    log_memoria("RAM ao final do teste")

    if ERROS:
        log("ERRO", f"{len(ERROS)} item(ns) critico(s) falharam")
        return 1
    log("OK", "ambiente Python pronto")
    return 0


if __name__ == "__main__":
    sys.exit(main())
