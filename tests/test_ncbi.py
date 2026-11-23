"""Testes de bioflora.ncbi sem rede: o Entrez e substituido por dublês."""

from __future__ import annotations

import io
import json

import pytest

from bioflora import ncbi

GENBANK_MINIMO = "LOCUS       TESTE0001  4 bp  DNA  circular PLN 01-JAN-2026\n//\n"


@pytest.fixture(autouse=True)
def dados_temporarios(monkeypatch, tmp_path):
    monkeypatch.setenv("BIOINFO_DADOS", str(tmp_path))
    monkeypatch.setenv("NCBI_EMAIL", "teste@exemplo.org")
    return tmp_path


def test_caminho_cache_dentro_do_data_lake(dados_temporarios):
    caminho = ncbi.caminho_cache("KY404932.1")
    assert caminho == dados_temporarios / "brutos" / "ncbi" / "KY404932.1.gb"


def test_formato_invalido():
    with pytest.raises(ValueError):
        ncbi.caminho_cache("KY404932.1", rettype="xml")


@pytest.mark.parametrize(
    ("accession", "esperado"),
    [("KY404932.1", True), ("NC_026450.1", True), ("KY404932", False)],
)
def test_accession_versionado(accession, esperado):
    assert ncbi.accession_versionado(accession) is esperado


def test_cache_evita_rede(monkeypatch):
    destino = ncbi.caminho_cache("TESTE0001.1")
    destino.write_text(GENBANK_MINIMO, encoding="utf-8")

    def proibido(**_kwargs):
        raise AssertionError("nao deveria acessar a rede")

    monkeypatch.setattr("Bio.Entrez.efetch", proibido)
    assert ncbi.baixar_registro("TESTE0001.1") == destino


def test_download_grava_arquivo_e_proveniencia(monkeypatch):
    monkeypatch.setattr("Bio.Entrez.efetch", lambda **_k: io.StringIO(GENBANK_MINIMO))
    destino = ncbi.baixar_registro("TESTE0001.1")
    assert destino.read_text(encoding="utf-8") == GENBANK_MINIMO
    meta = ncbi.ler_proveniencia(destino)
    assert meta["accession"] == "TESTE0001.1"
    assert meta["versionado"] is True
    assert len(meta["sha256"]) == 64
    assert not destino.with_name(destino.name + ".parcial").exists()
    json.dumps(meta)  # serializavel


def test_resposta_de_erro_nao_vira_cache(monkeypatch):
    monkeypatch.setattr(
        "Bio.Entrez.efetch", lambda **_k: io.StringIO("Error: ID list is empty")
    )
    with pytest.raises(ValueError):
        ncbi.baixar_registro("TESTE0002.1")
    assert not ncbi.caminho_cache("TESTE0002.1").exists()


def test_sem_email_falha_cedo(monkeypatch):
    monkeypatch.delenv("NCBI_EMAIL")
    with pytest.raises(RuntimeError):
        ncbi.baixar_registro("TESTE0003.1")
