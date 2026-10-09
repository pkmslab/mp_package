#' Assemble the final one-row feature table
#'
#' Internal helper used by [extract_mp_features()]. Combines the outputs
#' of every extraction/calculation step into a single-row table matching
#' the training dataset's feature columns, ready for user review and
#' manual editing before being passed to `predict_mp()`.
#'
#' @param resolved Output of [resolve_compound()].
#' @param props Output of [extract_physical_properties()].
#' @param pka_result Output of [extract_pka()].
#' @param fu_result Output of [extract_fu()].
#' @param drug_type Output of [classify_drug_type()].
#' @param logd Output of [calculate_logd()].
#' @param fraction_neutral Output of [calculate_fraction_neutral()].
#' @param bddcs_class Output of [lookup_bddcs()].
#'
#' @return A one-row data.frame with columns: Drug, SMILES, CAS,
#'   ChEMBL_ID, MW, LogP, PSA, HBD, HBA, pKa_acid, pKa_base, Fu, Type,
#'   LogD74, Fraction_neutral, BDDCS_class. Column-level source detail
#'   (which came from ChEMBL vs. PubChem, full pKa/Fu site detail) is
#'   attached as attributes (`property_sources`, `pka_detail`,
#'   `fu_detail`) rather than columns, since they don't collapse to a
#'   single value per drug.
#' @noRd
build_feature_table <- function(resolved, props, pka_result, fu_result,
                                drug_type, logd, fraction_neutral, bddcs_class) {

  table <- data.frame(
    Drug             = resolved$matched_name,
    SMILES           = resolved$smiles,
    CAS              = resolved$cas,
    ChEMBL_ID        = resolved$chembl_id,
    MW               = props$mw,
    LogP             = props$logp,
    PSA              = props$psa,
    HBD              = props$hbd,
    HBA              = props$hba,
    pKa_acid         = pka_result$pka_acid,
    pKa_base         = pka_result$pka_base,
    Fu               = fu_result$fu_median,
    Type             = drug_type,
    LogD74           = logd,
    Fraction_neutral = fraction_neutral,
    BDDCS_class      = bddcs_class,
    stringsAsFactors = FALSE
  )

  attr(table, "property_sources") <- props$sources
  attr(table, "pka_detail") <- pka_result$all_sites
  attr(table, "fu_detail") <- list(
    values  = fu_result$fu_values,
    n       = fu_result$n_values,
    source  = fu_result$source,
    flagged = fu_result$flagged
  )

  table
}
