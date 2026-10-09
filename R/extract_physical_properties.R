#' Extract physical properties for a compound
#'
#' Given a resolved compound (from [resolve_compound()]), extracts
#' molecular weight, LogP, polar surface area, and hydrogen bond donor/
#' acceptor counts. ChEMBL is used as the primary source when the compound
#' was resolved via ChEMBL; any individual property missing from ChEMBL
#' (or all properties, if the compound was resolved via PubChem) falls
#' back to PubChem, queried directly by SMILES structure.
#'
#' @param resolved The list returned by [resolve_compound()].
#'
#' @return A list with components:
#' \describe{
#'   \item{mw}{Molecular weight, or `NA`.}
#'   \item{logp}{Calculated LogP, or `NA`.}
#'   \item{psa}{Polar surface area, or `NA`.}
#'   \item{hbd}{Hydrogen bond donor count, or `NA`.}
#'   \item{hba}{Hydrogen bond acceptor count, or `NA`.}
#'   \item{sources}{A named list indicating which source (`"chembl"` or
#'     `"pubchem"`) each property came from, or `NA` if neither source had
#'     it.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' resolved <- resolve_compound("caffeine")
#' extract_physical_properties(resolved)
#' }
extract_physical_properties <- function(resolved) {

  if (missing(resolved) || !is.list(resolved) || is.null(resolved$smiles)) {
    stop("`resolved` must be the list returned by resolve_compound().", call. = FALSE)
  }

  chembl_props <- if (!is.null(resolved$chembl_id) && !is.na(resolved$chembl_id)) {
    query_chembl_properties(resolved$chembl_id)
  } else {
    empty_properties()
  }

  needed <- vapply(chembl_props[c("mw", "logp", "psa", "hbd", "hba")], is.na, logical(1))

  if (any(needed)) {
    message(
      if (all(needed)) "No physical properties found in ChEMBL. " else "Some properties missing from ChEMBL. ",
      "Checking PubChem for: ", paste(names(needed)[needed], collapse = ", "), "..."
    )
    pubchem_props <- query_pubchem_properties_by_smiles(resolved$smiles)
  } else {
    pubchem_props <- empty_properties()
  }

  final <- list(
    mw   = if (!is.na(chembl_props$mw))   chembl_props$mw   else pubchem_props$mw,
    logp = if (!is.na(chembl_props$logp)) chembl_props$logp else pubchem_props$logp,
    psa  = if (!is.na(chembl_props$psa))  chembl_props$psa  else pubchem_props$psa,
    hbd  = if (!is.na(chembl_props$hbd))  chembl_props$hbd  else pubchem_props$hbd,
    hba  = if (!is.na(chembl_props$hba))  chembl_props$hba  else pubchem_props$hba
  )

  sources <- list(
    mw   = property_source(chembl_props$mw,   pubchem_props$mw),
    logp = property_source(chembl_props$logp, pubchem_props$logp),
    psa  = property_source(chembl_props$psa,  pubchem_props$psa),
    hbd  = property_source(chembl_props$hbd,  pubchem_props$hbd),
    hba  = property_source(chembl_props$hba,  pubchem_props$hba)
  )

  still_missing <- vapply(final, is.na, logical(1))
  if (any(still_missing)) {
    message(
      "Could not find the following properties in either ChEMBL or PubChem: ",
      paste(names(still_missing)[still_missing], collapse = ", "), ". Returned as NA."
    )
  }

  c(final, list(sources = sources))
}

#' Determine which source a final property value came from
#' @noRd
property_source <- function(chembl_value, pubchem_value) {
  if (!is.na(chembl_value)) "chembl"
  else if (!is.na(pubchem_value)) "pubchem"
  else NA_character_
}

#' A properties list with all fields NA
#' @noRd
empty_properties <- function() {
  list(mw = NA_real_, logp = NA_real_, psa = NA_real_, hbd = NA_integer_, hba = NA_integer_)
}
