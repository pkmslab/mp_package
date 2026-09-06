#' Call pkasolver + classification pipeline on a single SMILES string
#'
#' Internal helper used by [extract_pka()]. Activates the target conda
#' environment, sources the bundled `pka_classification.py` module, and runs
#' the full site classification pipeline: charge-comparison direction,
#' SMARTS functional group identification, range filtering (pKa 3-11), and
#' selection of the lowest acidic / highest basic value.
#'
#' @param smiles A SMILES string.
#' @param envname Conda environment name.
#' @return A list with `pka_acid`, `pka_base` (numeric or NA), `fg_acid`,
#'   `fg_base` (character or NA), and `all_sites` (a data.frame of every
#'   site considered, including rejected ones, for transparency).
#' @noRd
call_pkasolver <- function(smiles, envname = "pka_prediction") {

  reticulate::use_condaenv(envname, required = TRUE)

  classification_module <- system.file(
    "python", "pka_classification.py",
    package = "mp.predict"
  )

  if (!nzchar(classification_module)) {
    stop(
      "Could not locate pka_classification.py inside the installed package. ",
      "Reinstall the package or check inst/python/ is present.",
      call. = FALSE
    )
  }

  pka_mod <- reticulate::import_from_path(
    "pka_classification",
    path = dirname(classification_module)
  )

  result <- tryCatch(
    pka_mod$extract_relevant_pka(smiles),
    error = function(e) {
      warning("pkasolver classification failed on '", smiles, "': ",
              conditionMessage(e), call. = FALSE)
      NULL
    }
  )

  if (is.null(result) || !is.null(result$error)) {
    if (!is.null(result) && identical(result$error, "invalid_smiles")) {
      warning("RDKit could not parse SMILES: '", smiles, "'.", call. = FALSE)
    }
    return(list(
      pka_acid = NA_real_, pka_base = NA_real_,
      fg_acid = NA_character_, fg_base = NA_character_,
      all_sites = data.frame()
    ))
  }

  all_sites <- sites_to_dataframe(result$all_sites)

  list(
    pka_acid  = if (is.null(result$pka_acid)) NA_real_ else result$pka_acid,
    pka_base  = if (is.null(result$pka_base)) NA_real_ else result$pka_base,
    fg_acid   = if (is.null(result$fg_acid)) NA_character_ else result$fg_acid,
    fg_base   = if (is.null(result$fg_base)) NA_character_ else result$fg_base,
    all_sites = all_sites
  )
}

#' Convert the Python list of dicts of per-site classifications into a
#' data.frame
#'
#' @param sites A Python list of dicts (as returned via reticulate), each
#'   with pka, direction, functional_group, method, kept.
#' @return A data.frame, one row per site. Empty data.frame if `sites` is
#'   empty.
#' @noRd
sites_to_dataframe <- function(sites) {
  if (length(sites) == 0) return(data.frame())

  data.frame(
    pka = vapply(sites, function(s) s$pka, numeric(1)),
    direction = vapply(sites, function(s) s$direction, character(1)),
    functional_group = vapply(
      sites, function(s) if (is.null(s$functional_group)) NA_character_ else s$functional_group,
      character(1)
    ),
    method = vapply(sites, function(s) s$method, character(1)),
    kept = vapply(sites, function(s) isTRUE(s$kept), logical(1)),
    stringsAsFactors = FALSE
  )
}
