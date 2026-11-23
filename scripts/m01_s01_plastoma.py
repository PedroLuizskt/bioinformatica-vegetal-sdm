"""Modulo 1 / Sessao 1 - o plastoma da araucaria como dado.

Etapas: ESearch com historico -> EFetch com cache e proveniencia -> SeqIO ->
tabela de anotacoes (Parquet) -> validacao das CDS contra /translation (Parquet).

Uso, a partir de 01_estudos com o .venv ativo:
    python scripts/m01_s01_plastoma.py                          # KY404932.1
    python scripts/m01_s01_plastoma.py --accession NC_026450.1  # outro plastoma
    python scripts/m01_s01_plastoma.py --offline                # so o cache local
    python scripts/m01_s01_plastoma.py --forcar                 # rebaixa do NCBI

Saidas em BIOINFO_DADOS (03_dados):
    brutos/ncbi/<accession>.gb (+ .gb.json de proveniencia)
    processados/plastomas/<accession>_features.parquet
    processados/plastomas/<accession>_cds_validacao.parquet
"""

from __future__ import annotations

import argparse
import sys

# bioflora antes de qualquer pacote cientifico (teto de threads do BLAS)
import bioflora  # noqa: F401
from bioflora.config import dir_dados
from bioflora.ncbi import (
    accession_versionado,
    baixar_registro,
    buscar,
    caminho_cache,
    ler_proveniencia,
)
from bioflora.recursos import log, log_memoria
from bioflora.sequencias import (
    features_para_dataframe,
    genes_repetidos,
    resumo_registro,
    validar_cds,
)

TERMO_PADRAO = (
    '"Araucaria angustifolia"[Organism] AND chloroplast[Title] '
    'AND "complete genome"[Title]'
)
MAX_SUSPEITAS = 20


def demonstrar_busca() -> None:
    """ESearch com servidor de historico; falha de rede nao interrompe o script."""
    try:
        r = buscar(TERMO_PADRAO)
    except Exception as exc:  # rede/NCBI fora: segue com o cache, se houver
        log("AVISO", f"ESearch indisponivel ({exc}); seguindo para o cache")
        return
    log("OK", f"ESearch: {r['total']} registro(s) | UIDs: {', '.join(r['ids'])}")
    log("OK", f"historico: query_key={r['query_key']} | WebEnv={r['webenv'][:24]}...")


def obter_arquivo(acc: str, offline: bool, forcar: bool):
    """Caminho do GenBank em cache (baixando se preciso) ou None em caso de erro."""
    if offline:
        arq = caminho_cache(acc)
        if not arq.exists():
            log("ERRO", f"{arq} nao existe; rode uma vez sem --offline")
            return None
        return arq
    ja_tinha = caminho_cache(acc).exists() and not forcar
    try:
        arq = baixar_registro(acc, forcar=forcar)
    except Exception as exc:
        log("ERRO", f"falha ao obter {acc}: {exc}")
        return None
    origem = "cache local" if ja_tinha else "NCBI (EFetch)"
    log("OK", f"{acc} via {origem}: {arq}")
    return arq


def relatar_validacao(val, arq_val) -> None:
    import pandas as pd  # importado depois de bioflora, de proposito

    n_validas = int(val["cds_valida"].sum())
    n_conferem = int(val["confere_anotacao"].fillna(False).sum())
    n_anotadas = int(val["tem_anotacao"].sum())
    log("OK", f"CDS: {len(val)} | validas pela tabela 11: {n_validas}")
    log("OK", f"traducao igual a /translation: {n_conferem} de {n_anotadas}")

    diverge = val["confere_anotacao"].eq(False).fillna(False)
    suspeitas = val[(~val["cds_valida"]) | diverge]
    for linha in suspeitas.head(MAX_SUSPEITAS).itertuples(index=False):
        partes = [linha.motivo] if linha.motivo else []
        if pd.notna(linha.diferencas_aa) and linha.diferencas_aa > 0:
            partes.append(f"{linha.diferencas_aa} aa diferentes da anotacao")
        detalhe = " | ".join(partes) or "verificar manualmente"
        log("AVISO", f"{linha.gene or '?':<8} inicio={linha.codon_inicio} | {detalhe}")
    if len(suspeitas) > MAX_SUSPEITAS:
        restantes = len(suspeitas) - MAX_SUSPEITAS
        log("AVISO", f"... e mais {restantes} (veja {arq_val.name})")


def main() -> int:
    parser = argparse.ArgumentParser(description="M1/S1 - plastoma como dado")
    parser.add_argument("--accession", default="KY404932.1")
    parser.add_argument("--offline", action="store_true", help="usa so o cache")
    parser.add_argument("--forcar", action="store_true", help="ignora o cache")
    args = parser.parse_args()

    acc = args.accession
    if not accession_versionado(acc):
        log("AVISO", f"{acc} sem versao (.1, .2...): o resultado pode mudar")

    # 1. Busca e download com cache
    if not args.offline:
        demonstrar_busca()
    arq = obter_arquivo(acc, args.offline, args.forcar)
    if arq is None:
        return 1
    prov = ler_proveniencia(arq)
    if prov:
        log("OK", f"baixado em {prov['baixado_em_utc']} | sha256 {prov['sha256'][:16]}")

    # 2. Leitura: um unico registro -> SeqIO.read
    from Bio import SeqIO

    record = SeqIO.read(arq, "genbank")
    r = resumo_registro(record)
    log("OK", f"{r['accession']} | {r['organismo']} | {r['descricao']}")
    log(
        "OK",
        f"{r['comprimento_pb']:,} pb | {r['topologia']} | GC {r['gc_pct']}% "
        f"| molecula {r['molecula']}",
    )
    log("OK", f"features por tipo: {r['n_features']}")
    log("OK", f"genes unicos: {r['genes_unicos']}")

    repetidos = genes_repetidos(record)
    if repetidos.empty:
        log("OK", "nenhum gene anotado em mais de uma copia")
    else:
        itens = repetidos.itertuples(index=False)
        lista = ", ".join(f"{g} ({t}, {n}x)" for t, g, n in itens)
        log("OK", f"genes com mais de uma copia: {lista}")

    # 3. Anotacoes -> Parquet
    saida = dir_dados() / "processados" / "plastomas"
    saida.mkdir(parents=True, exist_ok=True)
    feats = features_para_dataframe(record)
    arq_feats = saida / f"{acc}_features.parquet"
    feats.to_parquet(arq_feats, index=False)
    log("OK", f"{len(feats)} anotacoes -> {arq_feats}")

    # 4. Validacao das CDS contra a anotacao
    val = validar_cds(record)
    arq_val = saida / f"{acc}_cds_validacao.parquet"
    val.to_parquet(arq_val, index=False)
    relatar_validacao(val, arq_val)
    log("OK", f"validacao -> {arq_val}")

    log_memoria("RAM ao final")
    return 0


if __name__ == "__main__":
    sys.exit(main())
