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

import random

import numpy as np
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from bioflora.sequencias import (
    comparar_copias,
    coordenadas_genbank,
    features_para_dataframe,
    gc_janela,
    gc_janela_ingenua,
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

# Alfabeto IUPAC completo: a propriedade do exercicio 1 tambem vale para ambiguos.
IUPAC = "ACGTNRYSWKMBDHV"


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


def _registro_com_copias() -> SeqRecord:
    """Quatro copias de um 'tRNA' de 25 nt, uma para cada caso de comparar_copias.

    1: [0, 25) fita +       original
    2: [35, 60) fita -      complemento reverso do original (copia invertida)
    3: [70, 93) fita +      original sem as 2 primeiras bases (limite anotado)
    4: [103, 128) fita +    original com 1 substituicao interna
    """
    trna = "GCGGATTTAGCTCAGTTGGGAGAGC"
    mutado = trna[:10] + "A" + trna[11:]  # posicao 10 era C
    espaco = "T" * 10
    seq = (
        trna + espaco + str(Seq(trna).reverse_complement()) + espaco
        + trna[2:] + espaco + mutado
    )  # fmt: skip
    rec = SeqRecord(Seq(seq), id="SINT0002.1")
    locais = [(0, 25, 1), (35, 60, -1), (70, 93, 1), (103, 128, 1)]
    rec.features = [
        SeqFeature(FeatureLocation(i, f, strand=s), type="tRNA",
                   qualifiers={"gene": ["trnZ"]})
        for i, f, s in locais
    ]  # fmt: skip
    return rec


def test_comparar_copias_classifica_pares():
    pares = comparar_copias(_registro_com_copias()).set_index(["copia_a", "copia_b"])
    assert len(pares) == 6  # C(4, 2)

    p = pares.loc[(1, 2)]
    assert (p["relacao"], p["orientacao"]) == ("identica", "invertida")
    assert p["identidade_pct"] == 100.0

    p = pares.loc[(1, 3)]
    assert (p["relacao"], p["orientacao"]) == ("borda", "direta")
    assert (p["substituicoes"], p["gaps"]) == (0, 2)

    p = pares.loc[(1, 4)]
    assert p["relacao"] == "divergente"
    assert (p["substituicoes"], p["gaps"], p["identidades"]) == (1, 0, 24)


def test_comparar_copias_ignora_genes_de_copia_unica(registro):
    pares = comparar_copias(registro)
    assert set(pares["gene"]) == {"trnX"}


def test_resumo_registro(registro):
    r = resumo_registro(registro)
    assert r["comprimento_pb"] == 59
    assert r["genes_unicos"]["CDS"] == 5
    assert r["n_features"]["CDS"] == 5


def test_gc_percent_ignora_ambiguos():
    assert gc_percent("ATGC") == 50.0
    assert gc_percent("GGNN") == 100.0


def test_involucao_complemento_reverso_e_gc():
    """Exercicio 1 (Aula 1): rc(rc(s)) == s e GC(s) == GC(rc(s)).

    O complemento c troca G<->C e A<->T, entao leva o conjunto {G, C} nele mesmo;
    inverter a ordem nao muda contagens. Logo a contagem de G+C e a de bases
    validas sao as mesmas em s e em rc(s), e a igualdade e exata (mesmos inteiros
    na mesma divisao), sem tolerancia numerica. Vale tambem para os ambiguos:
    S (G/C) e W (A/T) sao fixos sob c, e R<->Y, K<->M, B<->V, D<->H sao
    descartados nos dois lados.
    """
    rng = random.Random(42)
    for _ in range(1000):
        s = "".join(rng.choices(IUPAC, k=rng.randint(10, 500)))
        rc = str(Seq(s).reverse_complement())
        assert str(Seq(rc).reverse_complement()) == s
        assert gc_percent(s) == gc_percent(rc)


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


@pytest.mark.parametrize(("janela", "passo"), [(50, 1), (50, 7), (200, 50), (30, 90)])
def test_gc_janela_igual_a_versao_ingenua(janela, passo):
    rng = np.random.default_rng(janela + passo)
    seq = "".join(rng.choice(list("ACGTNacgt"), size=2000))
    rapida = gc_janela(seq, janela=janela, passo=passo)
    ingenua = gc_janela_ingenua(seq, janela=janela, passo=passo)
    np.testing.assert_array_equal(rapida[0], ingenua[0])
    np.testing.assert_allclose(rapida[1], ingenua[1], equal_nan=True)


def test_gc_janela_ignora_caixa_e_aceita_seq():
    rng = np.random.default_rng(7)
    texto = "".join(rng.choice(list("ACGTN"), size=3000))
    referencia = gc_janela(texto, janela=300, passo=7)
    misturada = "".join(b.lower() if i % 3 == 0 else b for i, b in enumerate(texto))
    for variante in (misturada, Seq(texto), Seq(misturada)):
        centros, pct = gc_janela(variante, janela=300, passo=7)
        np.testing.assert_array_equal(centros, referencia[0])
        np.testing.assert_allclose(pct, referencia[1])


def test_gc_janela_sem_bases_validas_vira_nan():
    _, pct = gc_janela("NNNNGGGG", janela=4, passo=4)
    assert np.isnan(pct[0]) and pct[1] == 100.0


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
