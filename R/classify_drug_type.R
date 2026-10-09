#' Classify a drug's acid/base type from its extracted pKa values
#'
#' Determines whether a compound behaves as an acid, base, zwitterion, or
#' neutral species at physiological pH, based on which of `pka_acid` /
#' `pka_base` were found by [extract_pka()]. This classification feeds
#' into [calculate_logd()] and [calculate_fraction_neutral()], since the
#' LogD7.4 equation differs depending on drug type.
#'
#' @param pka_acid The acidic pKa value from [extract_pka()] (numeric or
#'   `NA`).
#' @param pka_base The basic pKa value from [extract_pka()] (numeric or
#'   `NA`).
#'
#' @return A single character string, one of:
#' \describe{
#'   \item{`"ACID"`}{`pka_acid` present, `pka_base` absent.}
#'   \item{`"BASE"`}{`pka_base` present, `pka_acid` absent.}
#'   \item{`"ZWITTERION"`}{Both `pka_acid` and `pka_base` present.}
#'   \item{`"NEUTRAL"`}{Neither present (no ionizable site found in the
#'     relevant range).}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- extract_pka("CC(=O)Oc1ccccc1C(=O)O")  # aspirin
#' classify_drug_type(result$pka_acid, result$pka_base)  # "ACID"
#' }
classify_drug_type <- function(pka_acid, pka_base) {

  has_acid <- !is.na(pka_acid)
  has_base <- !is.na(pka_base)

  if (has_acid && has_base) {
    "ZWITTERION"
  } else if (has_acid) {
    "ACID"
  } else if (has_base) {
    "BASE"
  } else {
    "NEUTRAL"
  }
}
