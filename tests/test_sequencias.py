"""Testes de bioflora.sequencias com um registro sintetico construido a mao.

Layout do registro (base 0, semiaberto):
    [0, 15)   CDS a  (+)  ATG AAA TTT GGG TAA            -> MKFG
    [15, 20)  espacador   CCCCC
    [20, 32)  CDS b  (-)  revcomp de ATG CCC AAA TAG     -> MPK
    [32, 38)  CDS c  (+)  exon 1: ATG GCT
    [38, 44)              intron: GTAAGT
    [44, 50)  CDS c  (+)  exon 2: GCT TAA                -> MAA
    [50, 59)  CDS e  (+)  ACG AAA TAA (inicio ACG, "editado" para M) -> anotado MK
    CDS d reaproveita a localizacao de a com anotacao errada (MKFA).
"""

from __future__ import annotations

import numpy as np
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from bioflora.sequencias import (
    coordenadas_genbank,
    features_para_dataframe,
    gc_janela,
    gc_percent,
    genes_repetidos,
    resumo_registro,
    validar_cds,
)

SEQ = (
    "ATGAAATTTGGGTAA"  # a
    "CCCCC"
    + str(Seq("ATGCCCAAATAG").reverse_complement())  # b (fita -)
    + "ATGGCT"
    + "GTAAGT"
    + "GCTTAA"  # c
    + "ACGAAATAA"  # e
)


def _cds(local, gene: str, traducao: str) -> SeqFeature:
    return SeqFeature(
        local, type="CDS", qualifiers={"gene": [gene], "translation": [traducao]}
    )


@pytest.fixture
def registro() -> SeqRecord:
    rec = SeqRecord(Seq(SEQ), id="SINT0001.1", name="SINT0001", description="sintetico")
    rec.annotations["molecule_type"] = "DNA"
    rec.annotations["topology"] = "circular"
    rec.annotations["organism"] = "Organismo sinteticus"
    rec.features = [
        _cds(FeatureLocation(0, 15, strand=1), "a", "MKFG"),
        _cds(FeatureLocation(20, 32, strand=-1), "b", "MPK"),
        _cds(
            CompoundLocation(
                [FeatureLocation(32, 38, strand=1), FeatureLocation(44, 50, strand=1)]
            ),
            "c",
            "MAA",
        ),
        _cds(FeatureLocation(0, 15, strand=1), "d", "MKFA"),
        _cds(FeatureLocation(50, 59, strand=1), "e", "MK"),
        SeqFeature(FeatureLocation(15, 20, strand=1), type="tRNA",
                   qualifiers={"gene": ["trnX"]}),
        SeqFeature(FeatureLocation(16, 20, strand=-1), type="tRNA",
                   qualifiers={"gene": ["trnX"]}),
    ]  # fmt: skip
    return rec


def test_comprimento_do_registro_sintetico():
    assert len(SEQ) == 59


def test_coordenadas_genbank_base1_fechado(registro):
    assert coordenadas_genbank(registro.features[0]) == (1, 15)
    assert coordenadas_genbank(registro.features[1]) == (21, 32)


def test_extract_respeita_fita_negativa(registro):
    assert str(registro.features[1].extract(registro.seq)) == "ATGCCCAAATAG"


def test_extract_remove_intron(registro):
    assert str(registro.features[2].extract(registro.seq)) == "ATGGCTGCTTAA"


def test_dataframe_de_anotacoes(registro):
    df = features_para_dataframe(registro)
    assert len(df) == 7
    c = df[df["gene"] == "c"].iloc[0]
    assert c["n_partes"] == 2
    assert c["comprimento"] == 12  # soma dos exons, sem o intron
    assert c["inicio"] == 33 and c["fim"] == 50
    assert df[df["gene"] == "b"].iloc[0]["fita"] == -1


def test_validar_cds(registro):
    val = validar_cds(registro).set_index("gene")
    for gene in ("a", "b", "c"):
        assert val.loc[gene, "cds_valida"]
        assert val.loc[gene, "confere_anotacao"]
        assert val.loc[gene, "diferencas_aa"] == 0
    assert not val.loc["d", "confere_anotacao"]
    assert val.loc["d", "diferencas_aa"] == 1
    assert not val.loc["e", "cds_valida"]
    assert val.loc["e", "codon_inicio"] == "ACG"
    assert val.loc["e", "diferencas_aa"] == 1  # assinatura de edicao de RNA


def test_genes_repetidos(registro):
    rep = genes_repetidos(registro)
    assert rep.to_dict("records") == [{"tipo": "tRNA", "gene": "trnX", "copias": 2}]


def test_resumo_registro(registro):
    r = resumo_registro(registro)
    assert r["comprimento_pb"] == 59
    assert r["genes_unicos"]["CDS"] == 5
    assert r["n_features"]["CDS"] == 5


def test_gc_percent_ignora_ambiguos():
    assert gc_percent("ATGC") == 50.0
    assert gc_percent("GGNN") == 100.0


def test_gc_janela_somas_de_prefixo():
    centros, pct = gc_janela("GGGGAAAA", janela=4, passo=4)
    np.testing.assert_allclose(centros, [2.5, 6.5])
    np.testing.assert_allclose(pct, [100.0, 0.0])


def test_gc_janela_igual_a_forca_bruta():
    rng = np.random.default_rng(42)
    seq = "".join(rng.choice(list("ACGTN"), size=5000))
    centros, pct = gc_janela(seq, janela=500, passo=250)
    for i, inicio in enumerate(range(0, len(seq) - 500 + 1, 250)):
        trecho = seq[inicio : inicio + 500]
        validas = sum(b in "ACGT" for b in trecho)
        esperado = 100 * sum(b in "GC" for b in trecho) / validas
        assert pct[i] == pytest.approx(esperado)
    assert len(centros) == len(pct)


def test_gc_janela_sequencia_curta():
    centros, pct = gc_janela("ACGT", janela=10)
    assert centros.size == 0 and pct.size == 0


def test_ida_e_volta_genbank(registro, tmp_path):
    arq = tmp_path / "sint.gb"
    SeqIO.write(registro, arq, "genbank")
    lido = SeqIO.read(arq, "genbank")
    assert str(lido.seq) == SEQ
    original = validar_cds(registro)
    relido = validar_cds(lido)
    assert original["cds_valida"].tolist() == relido["cds_valida"].tolist()
    assert features_para_dataframe(lido)["localizacao"].tolist() == (
        features_para_dataframe(registro)["localizacao"].tolist()
    )
