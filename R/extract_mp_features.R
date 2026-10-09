#' Extract all M/P prediction features for a drug
#'
#' The main entry point for feature extraction. Given a drug name or CAS
#' number, resolves the compound's structure, extracts physical properties,
#' pKa, unbound fraction (Fu), drug type, LogD, and BDDCS classification,
#' and returns a single-row, editable feature table -- the same columns
#' used to train the M/P prediction model. Review the returned table and
#' manually correct any value before passing it to `predict_mp()`.
#'
#' @param input A drug name (e.g. `"caffeine"`) or CAS number (e.g.
#'   `"58-08-2"`).
#' @param pH The pH to use for LogD/fraction-neutral calculations.
#'   Defaults to `7.4` (physiological plasma pH).
#'
#' @return A one-row data.frame -- see [build_feature_table()] for the
#'   full column list. Additional per-step detail (all pKa sites
#'   considered, all Fu values found, which source each physical property
#'   came from) is attached as attributes, viewable via `attributes()`.
#'
#' @details
#' Known gaps in the current implementation:
#' \itemize{
#'   \item Salt/hydrate stripping is not yet applied to compounds resolved
#'     via PubChem (ChEMBL's own search already tends to return parent
#'     structures in practice, but this has not been independently
#'     verified for every case).
#'   \item Fu extraction requires a ChEMBL ID; compounds resolved only via
#'     PubChem will have `Fu = NA`.
#'   \item pKa setup ([setup_pkasolver()]) must be run once per machine
#'     before this function can extract pKa.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' setup_pkasolver()  # once per machine
#' features <- extract_mp_features("caffeine")
#' features
#' }
extract_mp_features <- function(input, pH = 7.4) {

  message("Step 1/7: Resolving compound...")
  resolved <- resolve_compound(input)

  message("Step 2/7: Extracting physical properties...")
  props <- extract_physical_properties(resolved)

  message("Step 3/7: Extracting pKa...")
  pka_result <- extract_pka(resolved$smiles)

  message("Step 4/7: Extracting Fu...")
  if (!is.null(resolved$chembl_id) && !is.na(resolved$chembl_id)) {
    fu_result <- extract_fu(resolved$chembl_id)
  } else {
    message("No ChEMBL ID available for this compound (resolved via PubChem) -- Fu cannot be looked up. Returning NA.")
    fu_result <- list(
      fu_median = NA_real_, fu_values = numeric(0),
      n_values = 0L, source = NA_character_, flagged = FALSE
    )
  }

  message("Step 5/7: Classifying drug type and calculating LogD / fraction neutral...")
  drug_type <- classify_drug_type(pka_result$pka_acid, pka_result$pka_base)
  fraction_neutral <- calculate_fraction_neutral(pka_result, pH = pH)
  logd <- calculate_logd(props$logp, pka_result, pH = pH)

  message("Step 6/7: Looking up BDDCS class...")
  bddcs_class <- lookup_bddcs(resolved$matched_name)

  message("Step 7/7: Assembling feature table...")
  table <- build_feature_table(
    resolved = resolved, props = props, pka_result = pka_result,
    fu_result = fu_result, drug_type = drug_type, logd = logd,
    fraction_neutral = fraction_neutral, bddcs_class = bddcs_class
  )

  message("Done. Review the table below -- edit any value before using predict_mp().")
  table
}
