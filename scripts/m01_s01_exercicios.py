from pathlib import Path
import timeit
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from bioflora.sequencias import gc_janela

# Caminhos baseados na estrutura do projeto (01_estudos -> ../03_dados)
RAIZ_ESTUDOS = Path(__file__).resolve().parents[1]
DIR_DADOS = RAIZ_ESTUDOS.parent / "03_dados"
ARQ_GB_ANGUSTIFOLIA = DIR_DADOS / "brutos" / "ncbi" / "KY404932.1.gb"
DIR_OUT_PLASTOMAS = DIR_DADOS / "processados" / "plastomas"
DIR_OUT_PLASTOMAS.mkdir(parents=True, exist_ok=True)

record = SeqIO.read(ARQ_GB_ANGUSTIFOLIA, "genbank")

print("=" * 72)
print("EXERCICIO 2: Extracao de rbcL e matK (barcodes para o Modulo 3)")
print("=" * 72)
barcodes = []
for f in record.features:
    if f.type == "CDS":
        gene = f.qualifiers.get("gene", [""])[0]
        if gene in {"rbcL", "matK"}:
            seq_cds = f.extract(record.seq)
            n_pb = len(seq_cds)
            n_codons = n_pb // 3
            print(
                f"[OK] {gene}: {n_pb} pb | {n_codons} codons "
                f"(multiplo de 3: {n_pb % 3 == 0}) | inicio={seq_cds[:3]} | fim={seq_cds[-3:]}"
            )
            barcodes.append(
                SeqRecord(
                    seq_cds,
                    id=f"{record.id}_{gene}",
                    description=f"Araucaria angustifolia {gene} CDS ({n_pb} pb, {n_codons} codons)",
                )
            )

arq_fasta = DIR_OUT_PLASTOMAS / "KY404932.1_barcodes.fasta"
SeqIO.write(barcodes, arq_fasta, "fasta")
print(f"[OK] FASTA salvo em: {arq_fasta}\n")

print("=" * 72)
print("EXERCICIO 3: Verificacao de deposito (REFERENCE)")
print("=" * 72)
for i, ref in enumerate(record.annotations.get("references", []), start=1):
    print(f"--- Referencia {i} ---")
    print(f"Autores : {ref.authors}")
    print(f"Titulo  : {ref.title}")
    print(f"Journal : {ref.journal}")
    if ref.pubmed_id:
        print(f"PubMed  : {ref.pubmed_id}")
print()

print("=" * 72)
print("EXERCICIO 4: Investigacao dos genes com mais de uma copia (/gene vs /product)")
print("=" * 72)
alvos = {"rrn5", "trnD-GUC", "trnI-CAU", "trnR-UCU", "trnY-GUA"}
registros_rep = []
for f in record.features:
    if f.type in {"CDS", "tRNA", "rRNA"}:
        gene = f.qualifiers.get("gene", [""])[0]
        if gene in alvos:
            seq_ext = str(f.extract(record.seq))
            registros_rep.append(
                {
                    "tipo": f.type,
                    "gene": gene,
                    "product": f.qualifiers.get("product", ["-"])[0],
                    "note": f.qualifiers.get("note", ["-"])[0],
                    "localizacao": str(f.location),
                    "comprimento_pb": len(f.location),
                    "seq": seq_ext,
                }
            )

df_rep = pd.DataFrame(registros_rep)
for (tipo, gene), grupo in df_rep.groupby(["tipo", "gene"]):
    seqs_unicas = grupo["seq"].nunique()
    print(f"\nGene: {gene} ({tipo}) | Copias: {len(grupo)} | Sequencias distintas: {seqs_unicas}")
    for _, row in grupo.iterrows():
        print(
            f"  - loc={row['localizacao']} | {row['comprimento_pb']} pb | "
            f"/product='{row['product']}' | /note='{row['note']}'"
        )
print()

print("=" * 72)
print("EXERCICIO 5: Comparacao congenerica (KY404932.1 vs NC_026450.1)")
print("=" * 72)
arq_het_gb = DIR_DADOS / "brutos" / "ncbi" / "NC_026450.1.gb"
arq_het_feat = DIR_OUT_PLASTOMAS / "NC_026450.1_features.parquet"
arq_ang_feat = DIR_OUT_PLASTOMAS / "KY404932.1_features.parquet"

if arq_het_gb.exists() and arq_het_feat.exists():
    rec_het = SeqIO.read(arq_het_gb, "genbank")
    df_ang = pd.read_parquet(arq_ang_feat)
    df_het = pd.read_parquet(arq_het_feat)

    def resumo_comp(rec, df):
        sub = df[df["tipo"].isin(["CDS", "tRNA", "rRNA"])]
        cont = sub.groupby("tipo")["gene"].nunique().to_dict()
        gc = (rec.seq.count("G") + rec.seq.count("C")) / len(rec.seq) * 100
        return {
            "Especie": rec.annotations.get("organism", rec.id),
            "Accession": rec.id,
            "Comprimento (pb)": len(rec.seq),
            "GC (%)": round(gc, 2),
            "CDS unicas": cont.get("CDS", 0),
            "tRNA unicos": cont.get("tRNA", 0),
            "rRNA unicos": cont.get("rRNA", 0),
        }

    comp = pd.DataFrame([resumo_comp(record, df_ang), resumo_comp(rec_het, df_het)])
    print(comp.to_string(index=False))
else:
    print("[AVISO] Rode primeiro: python scripts\\m01_s01_plastoma.py --accession NC_026450.1")
print()

print("=" * 72)
print("EXERCICIO 6: Benchmark de complexidade gc_janela O(n) vs ingenua O(n*w)")
print("=" * 72)


def gc_janela_ingenua(seq, janela=1000, passo=250):
    s = str(seq).upper()
    n = len(s)
    centros = []
    valores = []
    for i in range(0, n - janela + 1, passo):
        fatia = s[i : i + janela]
        gc = (fatia.count("G") + fatia.count("C")) / janela * 100.0
        centros.append(i + janela / 2)
        valores.append(gc)
    return centros, valores


n_exec = 200
t_otimizado = timeit.timeit(
    lambda: gc_janela(record.seq, janela=1000, passo=250), number=n_exec
)
t_ingenuo = timeit.timeit(
    lambda: gc_janela_ingenua(record.seq, janela=1000, passo=250), number=n_exec
)

print(f"Plastoma: {len(record.seq):,} pb | janela w=1.000 | passo=250 | {n_exec} execucoes")
print(f"gc_janela (O(n))         : {t_otimizado:.4f} s ({t_otimizado/n_exec*1000:.3f} ms/exec)")
print(f"gc_janela_ingenua O(n*w) : {t_ingenuo:.4f} s ({t_ingenuo/n_exec*1000:.3f} ms/exec)")
print(f"Ganho de velocidade      : {t_ingenuo / t_otimizado:.1f}x mais rapida")