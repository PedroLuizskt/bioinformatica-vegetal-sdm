# -----------------------------------------------------------------------------
# Teste de fumaca do ambiente R (Aula 0).
#
# Uso, a partir de 01_estudos:
#   Rscript scripts/smoke_test.R            # inclui teste online (GBIF)
#   Rscript scripts/smoke_test.R --offline  # sem internet
#
# Sai com status 1 se algum item critico falhar.
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
offline <- "--offline" %in% args
erros <- 0L

log_msg <- function(nivel, ...) {
  cat(sprintf("[%s] %s\n", nivel, paste0(...)))
}
falha <- function(...) {
  log_msg("ERRO", ...)
  erros <<- erros + 1L
}

# 1. Sessao ---------------------------------------------------------------------------
# O CRAN so publica binarios Windows para a versao atual e a anterior do R; abaixo
# de 4.5, terra/sf/duckdb teriam de ser compilados do codigo-fonte.
if (getRversion() < "4.5.0") {
  falha(R.version.string, " e anterior a 4.5: rode o Passo 0 (instala o R atual) e o Passo 3 com -Reiniciar")
} else {
  log_msg("OK", R.version.string)
}
projeto <- tryCatch(renv::project(), error = function(e) NULL)
if (is.null(projeto)) {
  falha("renv nao esta ativo (rode o script a partir de 01_estudos)")
} else {
  log_msg("OK", "renv ativo, biblioteca: ", .libPaths()[1])
}

# 2. terra: gancho de limites + recorte + disco ---------------------------------------
tryCatch({
  suppressPackageStartupMessages(library(terra))
  log_msg("OK", "terra ", as.character(packageVersion("terra")),
          " | GDAL ", terra::gdal())

  opt <- terra::terraOptions(print = FALSE)
  tmp_esperado <- Sys.getenv("BIOINFO_TMP_TERRA", unset = "")
  if (is.list(opt) && !is.null(opt$memmax) && isTRUE(all.equal(as.numeric(opt$memmax), 4))) {
    log_msg("OK", "gancho do .Rprofile ativo: memmax = 4 GB, memfrac = ", opt$memfrac)
  } else {
    log_msg("AVISO", "nao confirmei memmax = 4 GB; confira o bloco no .Rprofile")
  }
  if (nzchar(tmp_esperado) && is.list(opt) && !is.null(opt$tempdir)) {
    mesmo <- normalizePath(opt$tempdir, winslash = "/", mustWork = FALSE) ==
      normalizePath(tmp_esperado, winslash = "/", mustWork = FALSE)
    if (isTRUE(mesmo)) {
      log_msg("OK", "tempdir do terra: ", opt$tempdir)
    } else {
      log_msg("AVISO", "tempdir do terra (", opt$tempdir, ") difere de BIOINFO_TMP_TERRA")
    }
  }

  # Raster global sintetico de 1 grau -> recorte para Minas Gerais -> disco -> leitura
  set.seed(42)
  global <- terra::rast(nrows = 180, ncols = 360, xmin = -180, xmax = 180,
                        ymin = -90, ymax = 90, crs = "EPSG:4326")
  terra::values(global) <- stats::runif(terra::ncell(global))
  ext_mg <- terra::ext(-51.1, -39.8, -22.95, -14.2)
  mg <- terra::crop(global, ext_mg)
  arq <- file.path(tempdir(), "smoke_mg.tif")
  terra::writeRaster(mg, arq, overwrite = TRUE)
  lido <- terra::rast(arq)
  log_msg("OK", "crop + writeRaster + rast: ", terra::ncell(lido), " celulas em MG (1 grau)")
  rm(global, mg, lido)
  terra::tmpFiles(remove = TRUE)
}, error = function(e) falha("terra: ", conditionMessage(e)))

# 3. sf + PROJ: SIRGAS 2000 geografico -> Policonica do Brasil ---------------------------
tryCatch({
  suppressPackageStartupMessages(library(sf))
  pt <- sf::st_sfc(sf::st_point(c(-45.28, -21.85)), crs = 4674)  # Cambuquira, MG
  proj <- sf::st_transform(pt, 5880)
  xy <- round(sf::st_coordinates(proj), 1)
  log_msg("OK", "sf ", as.character(packageVersion("sf")),
          " | EPSG:4674 -> EPSG:5880: (", xy[1], ", ", xy[2], ")")
}, error = function(e) falha("sf/PROJ: ", conditionMessage(e)))

# 4. DuckDB com teto de memoria + Parquet ----------------------------------------------
con <- NULL
tryCatch({
  # suppressMessages: o duckdb 2.x avisa sobre ~/.duckdb a cada conexao nova
  con <- suppressMessages(DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:"))
  DBI::dbExecute(con, "SET memory_limit = '2GB'")
  DBI::dbExecute(con, "SET threads = 4")
  pq <- gsub("\\\\", "/", file.path(tempdir(), "smoke.parquet"))
  df <- data.frame(especie = c("Araucaria angustifolia", "Handroanthus serratifolius"),
                   n = c(10L, 20L))
  duckdb::duckdb_register(con, "df", df)
  DBI::dbExecute(con, sprintf("COPY df TO '%s' (FORMAT parquet)", pq))
  total <- DBI::dbGetQuery(con, sprintf("SELECT sum(n) AS total FROM '%s'", pq))$total
  limite <- DBI::dbGetQuery(con, "SELECT current_setting('memory_limit') AS m")$m
  log_msg("OK", "duckdb ", as.character(packageVersion("duckdb")),
          " | Parquet ida e volta (soma = ", total, ") | memory_limit = ", limite)
}, error = function(e) falha("duckdb/Parquet: ", conditionMessage(e)),
finally = {
  if (!is.null(con)) try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
})

# 5. Pacotes de flora e do Quarto (knitr/rmarkdown) carregam -------------------------------
for (p in c("rgbif", "florabr", "CoordinateCleaner", "geodata", "knitr", "rmarkdown")) {
  ok <- requireNamespace(p, quietly = TRUE)
  if (ok) log_msg("OK", sprintf("%-18s %s", p, as.character(packageVersion(p))))
  else falha("pacote nao carrega: ", p)
}

# 6. reticulate -> Python do .venv (nao critico) ---------------------------------------------
tryCatch({
  py <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (!nzchar(py) || !file.exists(py)) {
    log_msg("AVISO", "RETICULATE_PYTHON nao aponta para um python.exe existente: '", py, "'")
  } else {
    bio <- reticulate::import("Bio")
    log_msg("OK", "reticulate -> ", py, " | Biopython ", bio$`__version__`)
  }
}, error = function(e) log_msg("AVISO", "reticulate: ", conditionMessage(e)))

# 7. GBIF online (nao critico) ----------------------------------------------------------------
if (offline) {
  log_msg("AVISO", "teste online do GBIF pulado (--offline)")
} else {
  tryCatch({
    nb <- rgbif::name_backbone(name = "Araucaria angustifolia", kingdom = "Plantae")
    chave <- nb$usageKey[1]
    n_br <- rgbif::occ_count(taxonKey = chave, country = "BR")
    log_msg("OK", "GBIF: Araucaria angustifolia usageKey = ", chave,
            " | ocorrencias no Brasil = ",
            formatC(as.numeric(n_br), format = "d", big.mark = ".", decimal.mark = ","))
  }, error = function(e) log_msg("AVISO", "GBIF indisponivel agora: ", conditionMessage(e)))
}

# 8. Memoria ---------------------------------------------------------------------------------
g <- gc()
log_msg("OK", sprintf("memoria R apos gc(): %.0f MB", sum(g[, 2])))

if (erros > 0L) {
  log_msg("ERRO", erros, " item(ns) critico(s) falharam")
  quit(save = "no", status = 1)
}
log_msg("OK", "ambiente R pronto")
