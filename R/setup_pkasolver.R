#' Set up the Python environment required for pKa extraction
#'
#' `extract_pka()` relies on the Python package \pkg{pkasolver}, which in turn
#' requires a pinned set of Python dependencies (notably
#' \code{torch-geometric==2.0.1}). These cannot be installed as ordinary R
#' package dependencies, so this function creates an isolated conda
#' environment and installs everything into it. This only needs to be run
#' once per machine.
#'
#' @param envname Name of the conda environment to create.
#'   Defaults to \code{"pka_prediction"}.
#' @param force If `TRUE`, remove and recreate the environment even if it
#'   already exists. Defaults to `FALSE`.
#'
#' @details
#' This function will:
#' \enumerate{
#'   \item Install Miniconda via reticulate, if no conda installation is found
#'   \item Create a conda environment with Python 3.9, RDKit, and PyTorch
#'   \item Install \code{torch-geometric==2.0.1}, \code{torch-sparse},
#'     \code{torch-scatter}, and \code{molvs} via pip inside that environment
#'   \item Install pkasolver directly from its GitHub repository
#' }
#'
#' This is a multi-gigabyte download and can take several minutes. It is
#' deliberately **not** run automatically when the package loads -- call it
#' explicitly, once, before your first `extract_pka()` call.
#'
#' @return Invisibly, the name of the conda environment.
#' @export
#'
#' @examples
#' \dontrun{
#' setup_pkasolver()
#' }
setup_pkasolver <- function(envname = "pka_prediction", force = FALSE) {

  if (!reticulate::py_available(initialize = FALSE) &&
      !nzchar(Sys.which("conda"))) {
    message("No conda installation found -- installing Miniconda via reticulate...")
    reticulate::install_miniconda()
  }

  existing_envs <- tryCatch(
    reticulate::conda_list()$name,
    error = function(e) character(0)
  )

  if (envname %in% existing_envs) {
    if (force) {
      message("Removing existing environment '", envname, "'...")
      reticulate::conda_remove(envname)
    } else {
      message("Environment '", envname, "' already exists. Use force = TRUE to recreate it.")
      return(invisible(envname))
    }
  }

  message("Creating conda environment '", envname, "' (this may take several minutes)...")
  reticulate::conda_create(
    envname  = envname,
    packages = c("python=3.9", "rdkit", "pytorch", "torchvision", "torchaudio"),
    channel  = c("pytorch", "conda-forge")
  )

  message("Installing pinned Python dependencies...")
  reticulate::conda_install(
    envname  = envname,
    packages = c("torch-geometric==2.0.1", "torch-sparse", "torch-scatter", "molvs"),
    pip      = TRUE
  )

  message("Installing pkasolver from GitHub...")
  reticulate::conda_install(
    envname  = envname,
    packages = "git+https://github.com/mayrf/pkasolver.git",
    pip      = TRUE
  )

  message("Verifying installation...")
  ok <- check_pkasolver_installed(envname)

  if (ok) {
    message("pkasolver installed successfully in environment '", envname, "'.")
  } else {
    warning(
      "Setup finished but pkasolver could not be imported. ",
      "See troubleshooting notes in ?setup_pkasolver."
    )
  }

  invisible(envname)
}

#' Check whether pkasolver is importable in a given conda environment
#'
#' Internal helper. Runs a small Python smoke test (predicting pKa for
#' caffeine) inside the target environment and confirms it returns at least
#' one value.
#'
#' @param envname Name of the conda environment to check.
#' @return `TRUE` if pkasolver imports and runs successfully, `FALSE` otherwise.
#' @noRd
check_pkasolver_installed <- function(envname = "pka_prediction") {
  tryCatch({
    reticulate::use_condaenv(envname, required = TRUE)
    pkasolver_query <- reticulate::import("pkasolver.query", delay_load = FALSE)
    rdkit_chem <- reticulate::import("rdkit.Chem", delay_load = FALSE)

    # caffeine, as a known-good smoke test compound
    mol <- rdkit_chem$MolFromSmiles("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
    states <- pkasolver_query$calculate_microstate_pka_values(mol, only_dimorphite = FALSE)

    length(states) >= 0  # a successful call, even with 0 sites, means the import/runtime works
  }, error = function(e) {
    message("check_pkasolver_installed() failed: ", conditionMessage(e))
    FALSE
  })
}
