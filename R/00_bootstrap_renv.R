# -----------------------------------------------------------------------------
# Aula 0 - bootstrap do renv no repositorio de estudos.
#
# Executar a partir de 01_estudos (o 03_setup_r.ps1 ja faz isso):
#   Rscript --no-init-file R/00_bootstrap_renv.R
#
# --no-init-file evita ler um .Rprofile antes de o renv existir.
# Idempotente: pode ser executado de novo sem efeitos colaterais.
# -----------------------------------------------------------------------------

log_msg <- function(nivel, ...) {
  cat(sprintf("[%s] %s\n", nivel, paste0(...)))
}

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  renv.consent = TRUE,
  Ncpus = 4L
)

if (!file.exists("_quarto.yml") || !dir.exists("R")) {
  log_msg("ERRO", "Execute este script a partir da pasta 01_estudos.")
  quit(save = "no", status = 1)
}

# 1. Biblioteca pessoal ---------------------------------------------------------
# Em modo nao interativo o R nao oferece criar a biblioteca pessoal sozinho;
# sem ela, install.packages() falha por falta de permissao em Program Files.
lib_user <- path.expand(Sys.getenv("R_LIBS_USER"))
if (!nzchar(lib_user)) {
  log_msg("ERRO", "R_LIBS_USER vazio; nao sei onde instalar o renv.")
  quit(save = "no", status = 1)
}
dir.create(lib_user, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib_user, .libPaths()))

# 2. renv na biblioteca pessoal ----------------------------------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  log_msg("OK", "Instalando renv em ", lib_user)
  install.packages("renv", lib = lib_user)
}
if (!requireNamespace("renv", quietly = TRUE)) {
  log_msg("ERRO", "Falha ao instalar o renv. Verifique a conexao com o CRAN.")
  quit(save = "no", status = 1)
}
log_msg("OK", "renv ", as.character(utils::packageVersion("renv")), " disponivel")

# 3. Inicializar o projeto (somente na primeira vez) --------------------------------
if (!file.exists(file.path("renv", "activate.R"))) {
  renv::init(
    project  = getwd(),
    bare     = TRUE,      # sem descobrir dependencias agora; instalamos no passo seguinte
    restart  = FALSE,
    settings = list(snapshot.type = "all")
  )
  log_msg("OK", "Projeto renv inicializado")
} else {
  log_msg("OK", "renv ja inicializado neste projeto")
}

# 4. Bloco de limites de recursos no .Rprofile --------------------------------------
# O renv cria o .Rprofile com source("renv/activate.R"). Acrescentamos o gancho que
# configura o terra em QUALQUER sessao deste projeto (console do VS Code, Rscript,
# Quarto). Detalhe importante: o terra recria suas opcoes no .onAttach, que roda
# DEPOIS do gancho "onLoad". Por isso registramos a mesma funcao em dois eventos:
#   - "attach": reaplica os limites depois de library(terra);
#   - "onLoad": cobre o uso so com terra::funcao(), sem library().
# Se o bloco ja existir, ele e substituido pela versao atual (atualizacao segura).
marcador_ini <- "# >>> bioinformatic: limites de recursos >>>"
marcador_fim <- "# <<< bioinformatic: limites de recursos <<<"
bloco <- c(
  marcador_ini,
  "# Gerado por R/00_bootstrap_renv.R. Teto de RAM e threads para o notebook de 8 GB.",
  "options(",
  "  repos = c(CRAN = \"https://cloud.r-project.org\"),",
  "  Ncpus = 4L,",
  "  mc.cores = 4L",
  ")",
  "local({",
  "  aplicar_limites_terra <- function(...) {",
  "    tmp <- Sys.getenv(\"BIOINFO_TMP_TERRA\", unset = \"\")",
  "    opcoes <- list(memfrac = 0.5, memmax = 4)",
  "    if (nzchar(tmp)) {",
  "      dir.create(tmp, recursive = TRUE, showWarnings = FALSE)",
  "      opcoes$tempdir <- tmp",
  "    }",
  "    do.call(terra::terraOptions, opcoes)",
  "  }",
  "  setHook(packageEvent(\"terra\", \"onLoad\"), aplicar_limites_terra)",
  "  setHook(packageEvent(\"terra\", \"attach\"), aplicar_limites_terra)",
  "})",
  marcador_fim
)

rprofile <- ".Rprofile"
atual <- if (file.exists(rprofile)) readLines(rprofile, warn = FALSE) else character()
ini <- which(atual == marcador_ini)
fim <- which(atual == marcador_fim)
if (length(ini) == 1 && length(fim) == 1 && fim > ini) {
  antes <- if (ini > 1) atual[seq_len(ini - 1)] else character()
  depois <- if (fim < length(atual)) atual[(fim + 1):length(atual)] else character()
  writeLines(c(antes, bloco, depois), rprofile)
  log_msg("OK", "Bloco de limites do .Rprofile atualizado")
} else {
  writeLines(c(atual, "", bloco), rprofile)
  log_msg("OK", "Bloco de limites adicionado ao .Rprofile")
}

log_msg("OK", "Bootstrap concluido")
