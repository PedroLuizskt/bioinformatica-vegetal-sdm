"""Guarda-corpos de recursos para o notebook de 8 GB (Aula 0).

Tres regras viram codigo aqui:
1. Teto de threads (4) para BLAS/OpenMP/numexpr.
2. Monitoramento de RAM do processo contra um teto seguro (~4,5 GB).
3. Conexao DuckDB ja nascida com memory_limit, threads e spill em disco.
"""

from __future__ import annotations

import os
import re
import sys
from typing import Any

MAX_THREADS = 4
TETO_RAM_GB = 4.5

# Lidas pelo BLAS/OpenMP no momento em que o numpy e carregado
VARIAVEIS_BLAS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
)
VARIAVEIS_THREADS = VARIAVEIS_BLAS + (
    "NUMEXPR_MAX_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)

_PADRAO_MEMORIA = re.compile(r"^\d+(\.\d+)?\s*(KB|MB|GB|KiB|MiB|GiB)$", re.IGNORECASE)


def log(nivel: str, mensagem: str) -> None:
    """Log em texto puro com marcadores [OK] / [AVISO] / [ERRO]."""
    print(f"[{nivel}] {mensagem}", flush=True)


def aplicar_limites_threads(
    n: int = MAX_THREADS, definidas_antes: set[str] | None = None
) -> list[str]:
    """Define o teto de threads nas variaveis de ambiente que ainda estiverem vazias.

    So tem efeito real se as variaveis BLAS existirem antes de o numpy ser
    carregado. ``definidas_antes`` e o conjunto de variaveis BLAS que ja estavam no
    ambiente antes de qualquer configuracao (ex.: vindas do .Renviron quando o
    Python roda pelo reticulate). Se o numpy ja estiver carregado e alguma delas
    faltava, o teto pode nao valer e um [AVISO] e emitido.
    """
    if definidas_antes is None:
        definidas_antes = {v for v in VARIAVEIS_BLAS if v in os.environ}
    n = max(1, min(int(n), MAX_THREADS))
    for var in VARIAVEIS_THREADS:
        os.environ.setdefault(var, str(n))
    os.environ.setdefault("GDAL_CACHEMAX", "512")

    avisos: list[str] = []
    faltavam = [v for v in VARIAVEIS_BLAS if v not in definidas_antes]
    ja_carregados = [m for m in ("numpy", "scipy") if m in sys.modules]
    if ja_carregados and faltavam:
        msg = (
            f"{', '.join(ja_carregados)} carregado(s) antes de bioflora sem "
            f"{', '.join(faltavam)} no ambiente; o teto de threads pode nao valer"
        )
        avisos.append(msg)
        log("AVISO", msg)
    return avisos


def memoria_status() -> dict[str, float]:
    """RAM do sistema e do processo atual, em GB."""
    import psutil

    vm = psutil.virtual_memory()
    rss = psutil.Process().memory_info().rss
    gib = 1024**3
    return {
        "total_gb": vm.total / gib,
        "disponivel_gb": vm.available / gib,
        "uso_sistema_pct": float(vm.percent),
        "processo_gb": rss / gib,
    }


def log_memoria(
    rotulo: str = "memoria", teto_gb: float = TETO_RAM_GB
) -> dict[str, float]:
    """Imprime e devolve o status de memoria; [AVISO] se o processo passar do teto."""
    s = memoria_status()
    nivel = "OK" if s["processo_gb"] <= teto_gb else "AVISO"
    log(
        nivel,
        f"{rotulo}: processo={s['processo_gb']:.2f} GB (teto {teto_gb:.1f}) | "
        f"disponivel={s['disponivel_gb']:.2f} de {s['total_gb']:.2f} GB",
    )
    return s


def conectar_duckdb(
    banco: str = ":memory:",
    memoria: str = "3GB",
    threads: int = MAX_THREADS,
    **config_extra: Any,
):
    """Abre uma conexao DuckDB com teto de memoria, threads e spill em disco.

    ``preserve_insertion_order = false`` reduz bastante a memoria em agregacoes
    e escritas grandes, ao custo de nao garantir a ordem original das linhas.
    """
    import duckdb

    from bioflora.config import dir_tmp

    if not _PADRAO_MEMORIA.match(memoria.strip()):
        raise ValueError(f"memoria invalida: {memoria!r} (use por exemplo '3GB')")

    tmp = dir_tmp() / "duckdb"
    tmp.mkdir(parents=True, exist_ok=True)
    config = {
        "memory_limit": memoria.strip(),
        "threads": max(1, min(int(threads), MAX_THREADS)),
        "temp_directory": tmp.as_posix(),
        "preserve_insertion_order": False,
    }
    config.update(config_extra)
    return duckdb.connect(banco, config=config)
