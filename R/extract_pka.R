#' Extract pKa values for a compound using pkasolver
#'
#' Predicts microstate pKa values for a molecule using the pkasolver Python
#' package (Mayr et al.), via a graph neural network model. A single
#' compound can have multiple ionizable sites, so this function may return
#' more than one pKa value.
#'
#' @param smiles A SMILES string for the (parent, salt-stripped) compound.
#' @param envname Name of the conda environment created by
#'   [setup_pkasolver()]. Defaults to `"pka_prediction"`.
#'
#' @return A list with components:
#' \describe{
#'   \item{pka_values}{Numeric vector of predicted pKa values (one per
#'     ionizable microstate). Empty if none are detected.}
#'   \item{n_sites}{Number of ionizable sites detected.}
#'   \item{pka_median}{Median of `pka_values`, or `NA` if none detected.}
#' }
#'
#' If more than one pKa value is found, a message points the user to
#' [select_pka_value()] to choose a specific value instead of the median.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' setup_pkasolver()  # once per machine
#' extract_pka("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")  # caffeine
#' }
extract_pka <- function(smiles, envname = "pka_prediction") {

  if (missing(smiles) || !is.character(smiles) || length(smiles) != 1 || !nzchar(smiles)) {
    stop("`smiles` must be a single non-empty character string.", call. = FALSE)
  }

  if (!check_pkasolver_installed(envname)) {
    stop(
      "pkasolver is not set up in conda environment '", envname, "'.\n",
      "Run setup_pkasolver() once before calling extract_pka().",
      call. = FALSE
    )
  }

  result <- call_pkasolver(smiles, envname = envname)

  if (length(result$pka_values) == 0) {
    message("No ionizable sites detected for this compound; pKa returned as NA.")
    return(list(pka_values = numeric(0), n_sites = 0L, pka_median = NA_real_))
  }

  if (length(result$pka_values) > 1) {
    message(
      length(result$pka_values), " pKa values found for this compound. ",
      "Returning the median. Use select_pka_value() to choose a specific value manually."
    )
  }

  list(
    pka_values = result$pka_values,
    n_sites    = result$n_sites,
    pka_median = stats::median(result$pka_values)
  )
}

#' Manually select a pKa value from multiple candidates
#'
#' Companion helper to [extract_pka()]. When a compound has multiple
#' predicted ionizable sites, use this to pick a specific pKa value (e.g. by
#' index or nearest to a known experimental value) instead of relying on the
#' default median.
#'
#' @param pka_result The list returned by [extract_pka()].
#' @param index Optional integer index into `pka_result$pka_values` to select
#'   directly.
#' @param value Optional numeric value to select if it exactly matches one of
#'   `pka_result$pka_values` (useful for scripting against known values).
#'
#' @return A single numeric pKa value.
#' @export
select_pka_value <- function(pka_result, index = NULL, value = NULL) {

  if (is.null(index) && is.null(value)) {
    stop("Provide either `index` or `value`.", call. = FALSE)
  }

  candidates <- pka_result$pka_values

  if (!is.null(index)) {
    if (index < 1 || index > length(candidates)) {
      stop("`index` is out of range for the available pKa values.", call. = FALSE)
    }
    return(candidates[index])
  }

  match_idx <- which(candidates == value)
  if (length(match_idx) == 0) {
    stop("`value` does not match any of the predicted pKa values.", call. = FALSE)
  }

  candidates[match_idx[1]]
}
