"""Acesso ao NCBI Entrez com cache local e proveniencia (Modulo 1).

O ``Bio.Entrez`` ja respeita o limite de taxa do NCBI (3 requisicoes/s sem chave,
10 com chave) e repete automaticamente falhas transitorias. Este modulo acrescenta
o que ele nao faz e que a reprodutibilidade exige:

- cache em disco por ``accession.version``: baixar uma vez, reproduzir sempre;
- arquivo de proveniencia ``.json`` ao lado de cada download (quando, o que, hash);
- validacao minima do conteudo, porque o NCBI pode devolver texto de erro com HTTP 200;
- escrita atomica (arquivo parcial + rename) para nunca deixar um cache corrompido.
"""

from __future__ import annotations

import hashlib
import json
import re
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from bioflora.config import configurar_entrez, dir_dados

# rettype -> (extensao do arquivo, prefixo esperado no conteudo)
FORMATOS: dict[str, tuple[str, str]] = {
    "gb": (".gb", "LOCUS"),
    "fasta": (".fasta", ">"),
}

_ACCESSION_VERSIONADO = re.compile(r"^[A-Z]{1,6}_?\d+\.\d+$")


def pasta_ncbi() -> Path:
    """Pasta de downloads brutos do NCBI dentro do data lake (criada se preciso)."""
    pasta = dir_dados() / "brutos" / "ncbi"
    pasta.mkdir(parents=True, exist_ok=True)
    return pasta


def caminho_cache(accession: str, rettype: str = "gb") -> Path:
    """Caminho do arquivo em cache para um accession e formato."""
    if rettype not in FORMATOS:
        raise ValueError(f"rettype nao suportado: {rettype!r} (use {sorted(FORMATOS)})")
    extensao, _ = FORMATOS[rettype]
    return pasta_ncbi() / f"{accession}{extensao}"


def accession_versionado(accession: str) -> bool:
    """True para identificadores com versao explicita (ex.: KY404932.1).

    Um accession sem versao aponta para a revisao mais recente do registro, que
    pode mudar; com versao, o conteudo e imutavel e o resultado e reproduzivel.
    """
    return bool(_ACCESSION_VERSIONADO.match(accession))


def _exigir_entrez() -> None:
    if not configurar_entrez():
        raise RuntimeError(
            "NCBI_EMAIL nao definido no .env; a politica do NCBI exige identificacao."
        )


def buscar(termo: str, db: str = "nucleotide", retmax: int = 20) -> dict[str, Any]:
    """ESearch com servidor de historico (usehistory='y').

    Devolve o total de registros, os primeiros ``retmax`` UIDs e o par
    ``WebEnv``/``query_key``, que permite buscar o conjunto inteiro em lotes
    depois, sem reenviar listas enormes de identificadores.
    """
    from Bio import Entrez

    _exigir_entrez()
    with Entrez.esearch(db=db, term=termo, retmax=retmax, usehistory="y") as handle:
        res = Entrez.read(handle)
    return {
        "termo": termo,
        "db": db,
        "total": int(res["Count"]),
        "ids": [str(i) for i in res["IdList"]],
        "webenv": str(res["WebEnv"]),
        "query_key": str(res["QueryKey"]),
    }


def _sha256(caminho: Path) -> str:
    h = hashlib.sha256()
    with caminho.open("rb") as f:
        for bloco in iter(lambda: f.read(1 << 20), b""):
            h.update(bloco)
    return h.hexdigest()


def baixar_registro(
    accession: str,
    rettype: str = "gb",
    db: str = "nucleotide",
    forcar: bool = False,
) -> Path:
    """Baixa um registro do NCBI para o cache local e devolve o caminho.

    Se o arquivo ja existir e ``forcar`` for False, nenhuma requisicao e feita.
    """
    destino = caminho_cache(accession, rettype)
    if destino.exists() and not forcar:
        return destino

    from Bio import Entrez

    _exigir_entrez()
    _, assinatura = FORMATOS[rettype]
    with Entrez.efetch(db=db, id=accession, rettype=rettype, retmode="text") as handle:
        texto = handle.read()
    if isinstance(texto, bytes):
        texto = texto.decode("utf-8")
    if not texto.lstrip().startswith(assinatura):
        raise ValueError(
            f"Resposta inesperada do NCBI para {accession}: {texto.strip()[:160]!r}"
        )

    parcial = destino.with_name(destino.name + ".parcial")
    parcial.write_text(texto, encoding="utf-8")
    parcial.replace(destino)

    import Bio

    proveniencia = {
        "accession": accession,
        "versionado": accession_versionado(accession),
        "db": db,
        "rettype": rettype,
        "baixado_em_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "bytes": destino.stat().st_size,
        "sha256": _sha256(destino),
        "fonte": "NCBI E-utilities (efetch)",
        "biopython": Bio.__version__,
    }
    destino.with_name(destino.name + ".json").write_text(
        json.dumps(proveniencia, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return destino


def ler_proveniencia(caminho: Path) -> dict[str, Any]:
    """Le o .json de proveniencia gravado ao lado de um download ({} se ausente)."""
    meta = caminho.with_name(caminho.name + ".json")
    if not meta.exists():
        return {}
    return json.loads(meta.read_text(encoding="utf-8"))
