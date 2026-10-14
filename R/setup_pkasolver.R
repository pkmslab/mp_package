#' Set up the Python environment required for pKa extraction
#'
#' `extract_pka()` relies on the Python package \pkg{pkasolver}, which in turn
#' requires a carefully matched set of Python dependencies (notably
#' \code{torch-geometric==2.0.1}, and compiled PyTorch extensions). These
#' cannot be installed as ordinary R package dependencies, so this function
#' creates an isolated conda environment and installs everything into it.
#' This only needs to be run once per machine.
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
#'   \item Create a conda environment with Python 3.9 and install RDKit
#'     (from conda-forge), then verify it imports
#'   \item Install a pinned, CPU-only PyTorch (2.3.1) from the pytorch channel
#'   \item Install \code{numpy<2} and \code{cairosvg} (with its native
#'     Cairo libraries) from conda-forge
#'   \item Install \code{torch-scatter} and \code{torch-sparse} from
#'     \strong{prebuilt wheels only} (no C++ compiler needed)
#'   \item Install \code{torch-geometric==2.0.1}, \code{molvs} and
#'     \code{svgutils} via pip
#'   \item Install pkasolver directly from its GitHub repository (it declares
#'     no dependencies of its own, which is why all of the above are listed
#'     explicitly)
#' }
#'
#' This is a multi-gigabyte download and can take 10-20 minutes. It is
#' deliberately **not** run automatically when the package loads -- call it
#' explicitly, once, before your first `extract_pka()` call.
#'
#' If the final check fails, restart R and run
#' `mp.predict:::check_pkasolver_installed()` once more: Python cannot swap
#' package versions inside an R session that has already loaded them.
#'
#' @return Invisibly, the name of the conda environment.
#' @export
#'
#' @examples
#' \dontrun{
#' setup_pkasolver()
#' }
setup_pkasolver <- function(envname = "pka_prediction", force = FALSE) {

  # Versions are pinned so that every machine gets the same, known-good stack.
  # torch-scatter / torch-sparse have no PyPI wheels on Windows, so they come
  # from the PyTorch Geometric wheel index, which is keyed by the torch
  # minor version (2.3.x -> "torch-2.3.0").
  torch_version       <- "2.3.1"
  torchvision_version <- "0.18.1"
  torchaudio_version  <- "2.3.1"
  pyg_wheel_index     <- "https://data.pyg.org/whl/torch-2.3.0+cpu.html"

  conda <- tryCatch(reticulate::conda_binary(), error = function(e) NULL)
  if (is.null(conda)) {
    message("No conda installation found -- installing Miniconda via reticulate...")
    reticulate::install_miniconda()
    conda <- reticulate::conda_binary()
  }

  existing_envs <- tryCatch(
    reticulate::conda_list(conda = conda)$name,
    error = function(e) character(0)
  )

  if (envname %in% existing_envs) {
    if (force) {
      message("Removing existing environment '", envname, "'...")
      reticulate::conda_remove(envname, conda = conda)
    } else {
      message("Environment '", envname, "' already exists. Use force = TRUE to recreate it.")
      return(invisible(envname))
    }
  }

  message("[1/6] Creating conda environment '", envname, "'...")
  reticulate::conda_create(
    envname  = envname,
    packages = "python=3.9",
    conda    = conda
  )

  message("[2/6] Installing RDKit (this may take a few minutes)...")
  reticulate::conda_install(
    envname  = envname,
    packages = "rdkit",
    channel  = "conda-forge",
    conda    = conda
  )

  if (!verify_rdkit_installed(envname, conda)) {
    stop(
      "RDKit could not be installed or imported in the '", envname, "' ",
      "conda environment. This can happen due to conda channel conflicts ",
      "on some platforms. Try running manually:\n",
      "  reticulate::conda_install(envname = \"", envname, "\", packages = \"rdkit\", channel = \"conda-forge\")",
      call. = FALSE
    )
  }

  message("[3/6] Installing PyTorch ", torch_version, " (CPU build; this may take several minutes)...")
  reticulate::conda_install(
    envname  = envname,
    packages = c(
      paste0("pytorch==", torch_version),
      paste0("torchvision==", torchvision_version),
      paste0("torchaudio==", torchaudio_version),
      "cpuonly"
    ),
    channel  = "pytorch",
    conda    = conda
  )

  message("[4/6] Installing numpy<2 and cairosvg...")
  reticulate::conda_install(
    envname  = envname,
    packages = c("numpy<2", "cairosvg"),
    channel  = "conda-forge",
    conda    = conda
  )

  message("[5/6] Installing torch-scatter and torch-sparse (prebuilt wheels)...")
  tryCatch(
    reticulate::conda_install(
      envname     = envname,
      packages    = c("torch-scatter", "torch-sparse"),
      pip         = TRUE,
      pip_options = c("--only-binary=:all:", "-f", pyg_wheel_index),
      conda       = conda
    ),
    error = function(e) {
      stop(
        "Could not find prebuilt torch-scatter / torch-sparse wheels for ",
        "this platform at ", pyg_wheel_index, ".\n",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  message("[6/6] Installing torch-geometric, molvs, svgutils and pkasolver...")
  reticulate::conda_install(
    envname  = envname,
    packages = c("torch-geometric==2.0.1", "molvs", "svgutils"),
    pip      = TRUE,
    conda    = conda
  )
  reticulate::conda_install(
    envname  = envname,
    packages = "git+https://github.com/mayrf/pkasolver.git",
    pip      = TRUE,
    conda    = conda
  )

  message("Verifying installation...")
  ok <- check_pkasolver_installed(envname)

  if (ok) {
    message("pkasolver installed successfully in environment '", envname, "'.")
  } else {
    warning(
      "Setup finished but the pkasolver check failed. Restart R and run ",
      "mp.predict:::check_pkasolver_installed() again; if it still fails, ",
      "see reticulate::py_last_error().",
      call. = FALSE
    )
  }

  invisible(envname)
}

#' Check whether RDKit specifically is importable in a conda environment
#'
#' Internal helper used by [setup_pkasolver()] to catch RDKit installation
#' failures early and with a clear message, rather than letting them
#' surface later as a cryptic Dimorphite-DL error at prediction time. This
#' was added after a cross-platform test found that installing RDKit
#' together with PyTorch in a single mixed-channel (`pytorch` +
#' `conda-forge`) conda_create() call could silently fail to resolve RDKit
#' on some platforms, with no error at install time.
#'
#' The check runs in a separate Python process rather than inside the R
#' session. Importing packages in-session during setup would load them into
#' memory before later install steps change their versions, and Python can't
#' reload them, leading to confusing "binary incompatibility" errors.
#'
#' @param envname Name of the conda environment to check.
#' @param conda Path to the conda binary used to create the environment.
#' @return `TRUE` if RDKit imports successfully, `FALSE` otherwise.
#' @noRd
verify_rdkit_installed <- function(envname, conda = "auto") {
  tryCatch({
    py <- reticulate::conda_python(envname, conda = conda)
    status <- suppressWarnings(
      system2(py, c("-c", shQuote("import rdkit.Chem")),
              stdout = FALSE, stderr = FALSE)
    )
    identical(as.integer(status), 0L)
  }, error = function(e) {
    FALSE
  })
}

#' Check whether pkasolver is importable and runs in a given conda environment
#'
#' Internal helper. Runs a small end-to-end smoke test (classifying the pKa
#' sites of caffeine through the package's own Python module) inside the
#' target environment and confirms it returns at least one site. Going
#' through the package's module, rather than calling pkasolver directly,
#' exercises the same Windows-safe subprocess handling that [extract_pka()]
#' uses.
#'
#' @param envname Name of the conda environment to check.
#' @return `TRUE` if pkasolver imports and runs successfully, `FALSE` otherwise.
#' @noRd
check_pkasolver_installed <- function(envname = "pka_prediction") {
  tryCatch({
    reticulate::use_condaenv(envname, required = TRUE)

    classification_module <- system.file(
      "python", "pka_classification.py", package = "mp.predict"
    )
    pka_mod <- reticulate::import_from_path(
      "pka_classification",
      path = dirname(classification_module),
      delay_load = FALSE
    )

    result <- NULL
    reticulate::py_capture_output({
      result <- pka_mod$extract_relevant_pka("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")
    })

    is.null(result$error) && length(result$all_sites) > 0
  }, error = function(e) {
    message("check_pkasolver_installed() failed: ", conditionMessage(e))
    FALSE
  })
}
