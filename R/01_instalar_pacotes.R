# -----------------------------------------------------------------------------
# Aula 0 - pacotes R essenciais e leves para comecar o Modulo 1.
#
# Executar a partir de 01_estudos (o 03_setup_r.ps1 ja faz isso):
#   Rscript R/01_instalar_pacotes.R
#
# O .Rprofile ativa o renv, entao tudo vai para a biblioteca do projeto.
#
# Desenho tolerante a falhas:
#   1. pre-checagem de acesso ao CRAN (falha cedo e com mensagem clara);
#   2. instala UM pacote por vez, pulando os ja instalados, para que um pacote
#      problematico nao impeca os outros;
#   3. no Windows prefere binarios (sem Rtools, nada compila);
#   4. registra tudo em 05_scratch/logs/r_pacotes_<data>.log;
#   5. faz renv::snapshot() do que instalou e sai com status 1 se algo faltou.
#
# Pacotes de SDM (maxnet, ENMeval, blockCV, biomod2, ecospat, usdm, GeoThinneR)
# entram nos modulos em que sao usados, para manter o ambiente enxuto.
# -----------------------------------------------------------------------------

raiz_ws <- Sys.getenv("BIOINFO_RAIZ", unset = normalizePath("..", winslash = "/"))
dir_logs <- file.path(raiz_ws, "05_scratch", "logs")
dir.create(dir_logs, recursive = TRUE, showWarnings = FALSE)
arq_log <- file.path(dir_logs, sprintf("r_pacotes_%s.log", format(Sys.time(), "%Y%m%d_%H%M%S")))

log_msg <- function(nivel, ...) {
  linha <- sprintf("[%s] %s", nivel, paste0(...))
  cat(linha, "\n", sep = "")
  cat(format(Sys.time(), "%H:%M:%S"), linha, "\n", file = arq_log, append = TRUE)
}

projeto <- tryCatch(renv::project(), error = function(e) NULL)
if (is.null(projeto)) {
  log_msg("ERRO", "renv nao esta ativo. Rode antes: Rscript --no-init-file R/00_bootstrap_renv.R")
  quit(save = "no", status = 1)
}
log_msg("OK", "Projeto renv ativo: ", projeto)
log_msg("OK", R.version.string, " | log em ", arq_log)

# Ordem importa: knitr/rmarkdown primeiro (Quarto), depois o pesado.
pacotes <- c(
  # Documentos Quarto (engine knitr) e ponte com Python
  "knitr", "rmarkdown", "here", "jsonlite", "reticulate",
  # Dados tabulares e interoperabilidade com o lado Python (Parquet via DuckDB)
  "DBI", "duckdb", "dplyr", "tidyr", "readr", "ggplot2",
  # Espacial (terra = raster/vetor em disco; sf = simple features; geodata = WorldClim/CMIP6)
  "terra", "sf", "geodata",
  # Ocorrencias e taxonomia de flora
  "rgbif", "florabr", "CoordinateCleaner",
  # Integracao com o VS Code (language server e visualizador de graficos)
  "languageserver", "httpgd"
)

no_windows <- identical(.Platform$OS.type, "windows")
tipo <- if (no_windows) "binary" else getOption("pkgType")
log_msg("OK", "Repositorio: ", paste(getOption("repos"), collapse = ", "), " | tipo preferido: ", tipo)

# 1. Pre-checagem de acesso ao CRAN -----------------------------------------------------
disponiveis <- tryCatch(
  utils::available.packages(type = tipo),
  error = function(e) {
    log_msg("ERRO", "Falha ao consultar o CRAN: ", conditionMessage(e))
    NULL
  }
)
if (is.null(disponiveis) || nrow(disponiveis) == 0) {
  log_msg("ERRO", "Sem acesso ao indice do CRAN. Verifique internet, proxy, VPN ou antivirus bloqueando o R.")
  quit(save = "no", status = 1)
}
log_msg("OK", nrow(disponiveis), " pacotes no indice ", tipo, " do CRAN")
sem_binario <- setdiff(pacotes, rownames(disponiveis))
if (length(sem_binario) > 0) {
  log_msg("AVISO", "Sem ", tipo, " no CRAN agora (tentara codigo-fonte): ", paste(sem_binario, collapse = ", "))
}

# 2. Instalacao pacote a pacote ---------------------------------------------------------------
falhas <- character()
for (p in pacotes) {
  if (requireNamespace(p, quietly = TRUE)) {
    log_msg("OK", sprintf("%-18s ja instalado (%s)", p, as.character(utils::packageVersion(p))))
    next
  }
  tipo_p <- if (p %in% sem_binario) "source" else tipo
  log_msg("OK", sprintf("%-18s instalando (%s) ...", p, tipo_p))
  erro <- tryCatch({
    renv::install(p, prompt = FALSE, type = tipo_p)
    NULL
  }, error = function(e) conditionMessage(e))

  if (is.null(erro) && requireNamespace(p, quietly = TRUE)) {
    log_msg("OK", sprintf("%-18s %s", p, as.character(utils::packageVersion(p))))
  } else {
    if (is.null(erro)) erro <- "instalou mas nao carrega"
    log_msg("ERRO", sprintf("%-18s %s", p, gsub("\\s+", " ", erro)))
    falhas <- c(falhas, p)
  }
  invisible(gc())
}

# 3. Trava do ambiente ----------------------------------------------------------------------
tryCatch({
  renv::snapshot(prompt = FALSE)
  log_msg("OK", "renv.lock atualizado")
}, error = function(e) log_msg("AVISO", "renv::snapshot falhou: ", conditionMessage(e)))

if (length(falhas) > 0) {
  log_msg("ERRO", length(falhas), " pacote(s) falharam: ", paste(falhas, collapse = ", "))
  log_msg("ERRO", "Rode este script de novo (ele pula o que ja instalou) e envie o log: ", arq_log)
  quit(save = "no", status = 1)
}
log_msg("OK", "Todos os ", length(pacotes), " pacotes essenciais prontos")
