"""Testes dos guarda-corpos de recursos e da configuracao do workspace."""

from __future__ import annotations

import os

import pytest

import bioflora
from bioflora import config, recursos


def test_versao_definida():
    assert bioflora.__version__ == "0.1.0"


def test_teto_de_threads_aplicado():
    for var in recursos.VARIAVEIS_THREADS:
        assert int(os.environ[var]) <= recursos.MAX_THREADS


def test_aplicar_limites_nao_sobrescreve(monkeypatch):
    monkeypatch.setenv("OMP_NUM_THREADS", "2")
    recursos.aplicar_limites_threads()
    assert os.environ["OMP_NUM_THREADS"] == "2"


def test_dir_dados_usa_variavel(monkeypatch, tmp_path):
    monkeypatch.setenv("BIOINFO_DADOS", str(tmp_path))
    assert config.dir_dados() == tmp_path


def test_dir_dados_fallback_para_repositorio(monkeypatch):
    monkeypatch.delenv("BIOINFO_DADOS", raising=False)
    assert config.dir_dados() == config.RAIZ / "dados"


def test_duckdb_nasce_com_limites(monkeypatch, tmp_path):
    monkeypatch.setenv("BIOINFO_TMP", str(tmp_path))
    con = recursos.conectar_duckdb(memoria="1GB", threads=8)
    try:
        threads = con.execute("SELECT current_setting('threads')").fetchone()[0]
        ordem = con.execute(
            "SELECT current_setting('preserve_insertion_order')"
        ).fetchone()[0]
    finally:
        con.close()
    assert int(threads) == recursos.MAX_THREADS  # 8 pedido, 4 aplicado
    assert ordem is False


def test_duckdb_rejeita_memoria_invalida():
    with pytest.raises(ValueError):
        recursos.conectar_duckdb(memoria="muito")


def test_memoria_status_coerente():
    s = recursos.memoria_status()
    assert 0 < s["processo_gb"] < s["total_gb"]


def test_cronometrar_mede_em_ms():
    import time

    rapido = recursos.cronometrar(lambda: sum(range(100)), repeticoes=2)
    lento = recursos.cronometrar(lambda: time.sleep(0.002), repeticoes=1)
    assert 0 < rapido < lento
    assert lento >= 2.0  # dormir 2 ms nunca custa menos de 2 ms


def test_cronometrar_rejeita_zero_repeticoes():
    with pytest.raises(ValueError):
        recursos.cronometrar(lambda: None, repeticoes=0)


def test_aviso_quando_numpy_veio_antes_sem_variaveis(monkeypatch):
    import numpy  # noqa: F401  (garante numpy em sys.modules)

    avisos = recursos.aplicar_limites_threads(definidas_antes=set())
    assert avisos and "OMP_NUM_THREADS" in avisos[0]


def test_sem_aviso_quando_variaveis_vieram_do_ambiente():
    import numpy  # noqa: F401

    avisos = recursos.aplicar_limites_threads(
        definidas_antes=set(recursos.VARIAVEIS_BLAS)
    )
    assert avisos == []
