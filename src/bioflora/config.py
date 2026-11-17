"""Caminhos do workspace e credenciais lidos do ambiente (.env).

Convencao do workspace:
- ``BIOINFO_DADOS`` aponta para o data lake local (03_dados). Sem a variavel
  (ex.: um clone limpo do repositorio), cai para ``<raiz do projeto>/dados``.
- ``BIOINFO_TMP`` aponta para temporarios descartaveis (05_scratch/tmp_py).
"""

from __future__ import annotations

import os
from pathlib import Path


def raiz_projeto() -> Path:
    """Sobe a partir deste arquivo ate achar o pyproject.toml do repositorio."""
    aqui = Path(__file__).resolve()
    for pai in aqui.parents:
        if (pai / "pyproject.toml").exists():
            return pai
    return Path.cwd()


RAIZ = raiz_projeto()


def carregar_env(arquivo: Path | None = None) -> bool:
    """Carrega o .env do projeto sem sobrescrever variaveis ja definidas.

    Retorna True se algum arquivo foi carregado.
    """
    try:
        from dotenv import load_dotenv
    except ImportError:
        return False
    alvo = arquivo if arquivo is not None else RAIZ / ".env"
    if not alvo.exists():
        return False
    return bool(load_dotenv(alvo, override=False))


def dir_dados() -> Path:
    """Diretorio do data lake (nao cria a pasta)."""
    valor = os.environ.get("BIOINFO_DADOS", "").strip()
    return Path(valor) if valor else RAIZ / "dados"


def dir_tmp() -> Path:
    """Diretorio de temporarios do Python (criado se nao existir)."""
    valor = os.environ.get("BIOINFO_TMP", "").strip()
    caminho = Path(valor) if valor else RAIZ / ".tmp"
    caminho.mkdir(parents=True, exist_ok=True)
    return caminho


def configurar_entrez(tool: str = "bioflora") -> bool:
    """Configura Bio.Entrez com email (obrigatorio) e chave de API (opcional).

    Retorna False se NCBI_EMAIL nao estiver definido.
    """
    from Bio import Entrez

    email = os.environ.get("NCBI_EMAIL", "").strip()
    if not email:
        return False
    Entrez.email = email
    Entrez.tool = tool
    chave = os.environ.get("NCBI_API_KEY", "").strip()
    if chave:
        Entrez.api_key = chave
    return True
