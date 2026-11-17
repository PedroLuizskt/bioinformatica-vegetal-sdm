# Bioinformática Vegetal e Modelagem de Distribuição de Espécies

[![CI](https://github.com/PedroLuizskt/bioinformatica-vegetal-sdm/actions/workflows/ci.yml/badge.svg)](https://github.com/PedroLuizskt/bioinformatica-vegetal-sdm/actions/workflows/ci.yml)
[![Licença: MIT](https://img.shields.io/badge/licen%C3%A7a-MIT-green.svg)](LICENSE)
![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)
![R 4.6](https://img.shields.io/badge/R-4.6-276DC3?logo=r&logoColor=white)
![Quarto](https://img.shields.io/badge/Quarto-apostila-75AADB?logo=quarto&logoColor=white)

Apostila aberta e laboratório reprodutível que unem **genômica comparativa vegetal
em Python** e **modelagem de distribuição de espécies (SDM/ENM) em R**, aplicadas a
árvores nativas do Brasil: ipês (*Handroanthus*), araucária (*Araucaria angustifolia*)
e jacarandás (*Dalbergia*).

A pergunta que orienta a trilha: **como ler a mesma espécie em três escalas —
sequência de DNA, linhagem evolutiva e paisagem — e usar isso para conservação?**

> Apostila renderizada: <https://pedroluizskt.github.io/bioinformatica-vegetal-sdm/>

---

## Roteiro

| Módulo | Tema | Python | R | Status |
|---|---|---|---|---|
| M0 | Setup e arquitetura | `venv`, pacote `bioflora`, Entrez | `renv`, limites de RAM no `terra` | concluído |
| M1 | Sequências e ocorrências como dados | `SeqIO`, plastomas, k-mers | GBIF, `florabr`, `CoordinateCleaner` | próximo |
| M2 | Alinhamento e espaço ambiental | Needleman–Wunsch, Smith–Waterman, BLAST+ | WorldClim, área M, colinearidade | planejado |
| M3 | Barcoding e algoritmos de nicho | MAFFT, barcode gap | `maxnet`, GAM, Random Forest | planejado |
| M4 | Filogenética e validação | Neighbor-Joining, Felsenstein, IQ-TREE 3 | `blockCV`, `ENMeval`, Boyce | planejado |
| M5 | Filogeografia e projeções | plastomas comparados, D de Tajima | `biomod2`, CMIP6, MESS | planejado |
| M6 | Filogenia × nicho | diversidade filogenética em H3 | sinal filogenético, `ecospat` | planejado |
| M7 | Capstone integrado | Ipês sob CMIP6: filogenia, nicho e perda de diversidade filogenética | | planejado |

Cada módulo gera uma nota de aula em Quarto (`modulos/`) e um incremento em um
projeto de portfólio independente.

## Estrutura

```text
.
|-- modulos/            notas de aula em Quarto, uma pasta por módulo
|-- src/bioflora/       pacote Python didático (guarda-corpos de recursos e,
|                       ao longo da trilha, algoritmos implementados do zero)
|-- tests/              testes pytest
|-- scripts/            testes de fumaça e utilitários de linha de comando
|-- R/                  bootstrap do renv e funções R reutilizáveis
|-- referencias/        bibliografia em BibTeX
|-- requirements/       dependências Python (intenção)
|-- requirements.lock.txt / renv.lock   travas exatas dos ambientes
`-- _quarto.yml         configuração do livro
```

## Princípios de engenharia

Todo o código roda em um notebook com **8 GB de RAM e sem GPU**:

- rasters sempre recortados para a área de estudo antes de qualquer extração;
- `terra` limitado a 4 GB e metade da RAM livre, com temporários em pasta descartável;
- no máximo 4 threads em BLAS/OpenMP e em qualquer paralelismo explícito;
- DuckDB com `memory_limit` e spill em disco para tabelas grandes;
- dados fora do Git, lidos da variável `BIOINFO_DADOS` (fallback: `./dados`);
- ambientes travados (`requirements.lock.txt`, `renv.lock`) e resultados congelados (`_freeze/`).

## Como reproduzir

Requisitos: Python 3.12, R ≥ 4.5 e Quarto.

```powershell
# Python (Windows; em Linux/macOS use requirements/base.txt no lugar do lock)
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.lock.txt
.\.venv\Scripts\python.exe -m pip install --no-deps -e .

# R: o .Rprofile ativa o renv; restore instala as versões travadas
Rscript -e "renv::restore(prompt = FALSE)"

# Verificar e renderizar
.\.venv\Scripts\python.exe scripts\smoke_test.py
.\.venv\Scripts\python.exe -m pytest
Rscript scripts\smoke_test.R
quarto render
```

Para consultas ao NCBI, crie um `.env` a partir de `.env.example` e preencha `NCBI_EMAIL`.

## Como citar

Os metadados de citação estão em [`CITATION.cff`](CITATION.cff); o GitHub mostra o
botão *Cite this repository* a partir dele.

## Autor

**Pedro Luiz Rodrigues Vaz de Melo**: Engenheiro Florestal e Cientista de Dados Geoespaciais, Cambuquira (MG).

- GitHub: [@PedroLuizskt](https://github.com/PedroLuizskt)
- Portfólio: <https://pedroluizskt.github.io/>
<!-- Para exibir, preencha e mova para fora deste comentário:
- LinkedIn: https://www.linkedin.com/in/SEU-PERFIL
- Email: seu.email@exemplo.com
-->

## Licença

Código sob [licença MIT](LICENSE). Textos das notas de aula sob
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/deed.pt-br).
