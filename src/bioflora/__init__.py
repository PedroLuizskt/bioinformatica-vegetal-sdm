"""bioflora: pacote didatico da trilha Bioinformatica Vegetal (Python) + SDM (R).

Regra de uso em scripts pesados: ``import bioflora`` ANTES de numpy/pandas.
O import carrega o ``.env`` do projeto e fixa o teto de threads de BLAS/OpenMP,
e essas variaveis so tem efeito se existirem antes de o numpy ser carregado.
"""

import os

from bioflora.config import carregar_env
from bioflora.recursos import VARIAVEIS_BLAS, aplicar_limites_threads

# 0) o que ja vinha do ambiente (ex.: .Renviron quando o Python roda pelo reticulate);
# 1) .env do projeto; 2) defaults para o que ainda faltar.
_blas_definidas_antes = {v for v in VARIAVEIS_BLAS if v in os.environ}
carregar_env()
aplicar_limites_threads(definidas_antes=_blas_definidas_antes)

__version__ = "0.1.0"

__all__ = ["__version__"]
