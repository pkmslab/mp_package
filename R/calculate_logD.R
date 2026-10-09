#' Calculate fraction neutral at a given pH
#'
#' Computes the fraction of a compound present in its neutral (un-ionized)
#' form at a given pH, using ALL valid ionizable sites found by
#' [extract_pka()] (not just the collapsed single acid/base values) --
#' correctly handling drugs with two acidic or two basic sites via the
#' standard polyprotic dissociation equations, matching the method used to
#' build the training dataset.
#'
#' @param pka_result The list returned by [extract_pka()] (must include
#'   `all_sites`).
#' @param pH The pH to calculate at. Defaults to `7.4` (physiological
#'   plasma pH).
#'
#' @return A single numeric value: the fraction neutral. `1` if no valid
#'   ionizable sites were found (compound treated as fully neutral).
#'
#' @details
#' Sites are taken from `pka_result$all_sites` where `kept == TRUE`, split
#' into acidic (sorted ascending) and basic (sorted descending) lists.
#' \itemize{
#'   \item 0 sites on a side: fraction for that side = 1
#'   \item 1 site: standard monoprotic equation
#'   \item 2+ sites: the two most relevant pKas (closest to being fully
#'     ionized/deprotonated first) are used in the diprotic equation; any
#'     additional sites beyond the first two are not currently used
#' }
#' Diprotic acid: \eqn{f_{acid} = 1 / (1 + 10^{(pH-pKa_1)} + 10^{(2pH-pKa_1-pKa_2)})}
#'
#' Diprotic base: \eqn{f_{base} = 1 / (1 + 10^{(pKa_1-pH)} + 10^{(pKa_1+pKa_2-2pH)})}
#'
#' Note: `all_sites` is filtered upstream in [extract_pka()] to pKa range
#' 3-11 with an approved functional group -- this is a different (and
#' stricter) filter than the `abs(pKa - pH) <= 4` filter used in the
#' original raw-pKasolver-based diprotic calculation. This may produce
#' slightly different results for edge-case pKas (e.g. between 3-3.4 or
#' 11-11.4) and should be confirmed against the training methodology.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pka_result <- extract_pka("CC(=O)Nc1ccc(O)cc1")  # acetaminophen
#' calculate_fraction_neutral(pka_result)
#' }
calculate_fraction_neutral <- function(pka_result, pH = 7.4) {

  if (is.null(pka_result$all_sites) || nrow(pka_result$all_sites) == 0) {
    return(1)
  }

  sites <- pka_result$all_sites
  kept_sites <- sites[isTRUE(sites$kept) | sites$kept == TRUE, , drop = FALSE]

  if (nrow(kept_sites) == 0) {
    return(1)
  }

  acid_pkas <- sort(kept_sites$pka[kept_sites$direction == "ACID"])
  base_pkas <- sort(kept_sites$pka[kept_sites$direction == "BASE"], decreasing = TRUE)

  f_acid <- fraction_neutral_side(acid_pkas, pH, is_acid = TRUE)
  f_base <- fraction_neutral_side(base_pkas, pH, is_acid = FALSE)

  f_acid * f_base
}

#' Compute the fraction-neutral contribution from one side (acid or base)
#'
#' Internal helper for [calculate_fraction_neutral()]. Applies the mono- or
#' di-protic equation depending on how many valid pKas are on this side.
#'
#' @param pkas A sorted numeric vector of pKa values for this side (acid:
#'   ascending: base: descending).
#' @param pH The pH to calculate at.
#' @param is_acid `TRUE` for the acidic side, `FALSE` for the basic side --
#'   determines which form of the equation to use.
#' @return A single numeric fraction (0-1).
#' @noRd
fraction_neutral_side <- function(pkas, pH, is_acid) {
  n <- length(pkas)

  if (n == 0) {
    return(1)
  }

  if (n == 1) {
    if (is_acid) {
      return(1 / (1 + 10^(pH - pkas[1])))
    } else {
      return(1 / (1 + 10^(pkas[1] - pH)))
    }
  }

  # 2 or more: use the two most relevant pKas (diprotic equation).
  # Any sites beyond the first two are not currently incorporated.
  pka1 <- pkas[1]
  pka2 <- pkas[2]

  if (is_acid) {
    1 / (1 + 10^(pH - pka1) + 10^(2 * pH - pka1 - pka2))
  } else {
    1 / (1 + 10^(pka1 - pH) + 10^(pka1 + pka2 - 2 * pH))
  }
}

#' Calculate LogD at a given pH
#'
#' Computes the distribution coefficient (LogD) from LogP and the fraction
#' neutral (via [calculate_fraction_neutral()], which handles diprotic
#' sites correctly) at a given pH.
#' \deqn{LogD = LogP + \log_{10}(f_{neutral})}
#'
#' @param logp Calculated LogP (e.g. from [extract_physical_properties()]).
#' @param pka_result The list returned by [extract_pka()].
#' @param pH The pH to calculate at. Defaults to `7.4`.
#'
#' @return A single numeric LogD value, or `NA` if `logp` is `NA`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pka_result <- extract_pka("CC(=O)Nc1ccc(O)cc1")  # acetaminophen
#' calculate_logd(logp = 0.46, pka_result = pka_result)
#' }
calculate_logd <- function(logp, pka_result, pH = 7.4) {

  if (is.na(logp)) {
    return(NA_real_)
  }

  f_neutral <- calculate_fraction_neutral(pka_result, pH)

  logp + log10(f_neutral)
}
