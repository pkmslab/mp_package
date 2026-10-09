#' Extract unbound fraction in plasma (Fu) for a compound
#'
#' Queries ChEMBL's activity data for human plasma Fu measurements, falling
#' back to plasma protein binding (PPB) data (Fu = 1 - PPB/100) if no direct
#' Fu measurement is available. Returns NA if neither is found
#'
#' @param chembl_id A ChEMBL molecule ID (e.g. from [resolve_compound()]).
#'
#' @return A list with components:
#' \describe{
#'   \item{fu_median}{The median Fu value across all valid measurements
#'     found, or `NA` if none.}
#'   \item{fu_values}{Numeric vector of all individual valid Fu values used
#'     to compute the median (already converted from PPB if that was the
#'     source).}
#'   \item{n_values}{Number of valid measurements found.}
#'   \item{source}{One of `"chembl_fu"`, `"chembl_ppb"`, or `NA` if nothing
#'     was found.}
#'   \item{flagged}{`TRUE` if any individual raw value was exactly 0 or 1
#'     (direct Fu measurements only).}
#' }
#'
#' If more than one valid value is found, a message points the user to
#' [select_fu_value()] to choose a specific value instead of the median.
#'
#' @details
#' Filtering applied to ChEMBL activity records (matches the method used
#' to build the training dataset):
#' \itemize{
#'   \item Human only: `target_organism` is `NA` or `"Homo sapiens"`
#'   \item Animal keywords excluded from `assay_description` (rat, mouse,
#'     mice, animal, rabbit, dog, monkey -- word-boundary matched)
#'   \item `"tissue"` excluded from `assay_description`
#'   \item `assay_description` must contain `"plasma"` or `"human"`
#'   \item Direct Fu values restricted to the range 0-1
#'   \item PPB values (fallback only) restricted to the range 0-100
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' extract_fu("CHEMBL25")  # aspirin
#' }
extract_fu <- function(chembl_id) {

  if (missing(chembl_id) || !is.character(chembl_id) || length(chembl_id) != 1 || !nzchar(chembl_id)) {
    stop("`chembl_id` must be a single non-empty character string.", call. = FALSE)
  }

  fu_result <- query_chembl_fu(chembl_id)

  if (fu_result$n_values > 0) {
    result <- list(
      fu_median = stats::median(fu_result$values),
      fu_values = fu_result$values,
      n_values  = fu_result$n_values,
      source    = "chembl_fu",
      flagged   = fu_result$flagged
    )
  } else {
    message("No direct Fu data found in ChEMBL. Checking plasma protein binding (PPB) data...")
    ppb_result <- query_chembl_ppb(chembl_id)

    if (ppb_result$n_values > 0) {
      result <- list(
        fu_median = stats::median(ppb_result$fu_values),
        fu_values = ppb_result$fu_values,
        n_values  = ppb_result$n_values,
        source    = "chembl_ppb",
        flagged   = FALSE  # the 0/1 flag only applies to direct Fu measurements
      )
    } else {
      message("No Fu or PPB data found in ChEMBL. Returning NA for manual entry.")
      result <- list(
        fu_median = NA_real_,
        fu_values = numeric(0),
        n_values  = 0L,
        source    = NA_character_,
        flagged   = FALSE
      )
    }
  }

  if (result$n_values > 1) {
    message(
      result$n_values, " valid Fu values found for this compound. ",
      "Returning the median. Use select_fu_value() to choose a specific value manually."
    )
  }

  if (isTRUE(result$flagged)) {
    message("Note: at least one raw Fu value was exactly 0 or 1 : flagged for manual literature review.")
  }

  result
}

#' Manually select an Fu value from multiple candidates
#'
#' Companion helper to [extract_fu()]. When a compound has multiple valid
#' Fu measurements, use this to pick a specific value (e.g. by index) instead
#' of relying on the default median.
#'
#' @param fu_result The list returned by [extract_fu()].
#' @param index Optional integer index into `fu_result$fu_values` to select
#'   directly.
#' @param value Optional numeric value to select if it exactly matches one of
#'   `fu_result$fu_values`.
#'
#' @return A single numeric Fu value.
#' @export
select_fu_value <- function(fu_result, index = NULL, value = NULL) {

  if (is.null(index) && is.null(value)) {
    stop("Provide either `index` or `value`.", call. = FALSE)
  }

  candidates <- fu_result$fu_values

  if (!is.null(index)) {
    if (index < 1 || index > length(candidates)) {
      stop("`index` is out of range for the available Fu values.", call. = FALSE)
    }
    return(candidates[index])
  }

  match_idx <- which(candidates == value)
  if (length(match_idx) == 0) {
    stop("`value` does not match any of the available Fu values.", call. = FALSE)
  }

  candidates[match_idx[1]]
}
