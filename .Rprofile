source("renv/activate.R")

# >>> bioinformatic: limites de recursos >>>
# Gerado por R/00_bootstrap_renv.R. Teto de RAM e threads para o notebook de 8 GB.
options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  Ncpus = 4L,
  mc.cores = 4L
)
local({
  aplicar_limites_terra <- function(...) {
    tmp <- Sys.getenv("BIOINFO_TMP_TERRA", unset = "")
    opcoes <- list(memfrac = 0.5, memmax = 4)
    if (nzchar(tmp)) {
      dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
      opcoes$tempdir <- tmp
    }
    do.call(terra::terraOptions, opcoes)
  }
  setHook(packageEvent("terra", "onLoad"), aplicar_limites_terra)
  setHook(packageEvent("terra", "attach"), aplicar_limites_terra)
})
# <<< bioinformatic: limites de recursos <<<
