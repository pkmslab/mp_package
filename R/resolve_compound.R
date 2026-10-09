#' Resolve a drug name or CAS number to a canonical structure
#'
#' Looks up a compound by name or CAS registry number, querying ChEMBL first
#' and falling back to PubChem if no match is found. Returns the parent
#' (salt-stripped) canonical SMILES along with identifiers from whichever
#' source matched.
#'
#' @param input A drug name (e.g. `"caffeine"`) or CAS number
#'   (e.g. `"58-08-2"`). CAS format is detected automatically.
#'
#' @return A list with components:
#' \describe{
#'   \item{smiles}{Canonical SMILES string of the parent structure.}
#'   \item{cas}{CAS registry number, if available (from either source).}
#'   \item{chembl_id}{ChEMBL ID, if resolved via ChEMBL (`NA` otherwise).}
#'   \item{pubchem_cid}{PubChem CID, if resolved via PubChem (`NA` otherwise).}
#'   \item{source}{Either `"chembl"` or `"pubchem"`, indicating which
#'     database the match came from.}
#'   \item{matched_name}{The name/identifier as recorded by the source
#'     database (useful for confirming the right compound was found).}
#' }
#'
#' @details
#' ChEMBL's parent molecule structure is used directly where available,
#' since ChEMBL already resolves salts/hydrates to a parent compound
#' internally. PubChem-sourced results are NOT guaranteed to be
#' salt-stripped and should be passed through [strip_salts_and_hydrates()]
#' downstream.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' resolve_compound("caffeine")
#' resolve_compound("58-08-2")
#' }
resolve_compound <- function(input) {

  if (missing(input) || !is.character(input) || length(input) != 1 || !nzchar(trimws(input))) {
    stop("`input` must be a single non-empty character string (drug name or CAS number).", call. = FALSE)
  }

  input <- trimws(input)
  is_cas <- is_cas_number(input)

  message(
    "Resolving '", input, "' as a ", if (is_cas) "CAS number" else "drug name",
    "  querying ChEMBL..."
  )

  result <- if (is_cas) query_chembl_by_cas(input) else query_chembl_by_name(input)

  if (!is.null(result)) {
    message("Match found in ChEMBL: ", result$matched_name)
    return(fill_cas_via_pubchem(result))
  }

  message("No match in ChEMBL. Falling back to PubChem...")

  result <- if (is_cas) query_pubchem_by_cas(input) else query_pubchem_by_name(input)

  if (!is.null(result)) {
    message("Match found in PubChem: ", result$matched_name)
    return(fill_cas_via_pubchem(result))
  }

  stop(
    "Could not resolve '", input, "' in either ChEMBL or PubChem. ",
    "Check spelling, or try the CAS number if you have one.",
    call. = FALSE
  )
}

#' Ensure a resolved compound has a CAS number, via PubChem if needed
#'
#' Internal helper. ChEMBL's molecule_synonyms was found to not reliably
#' curate CAS numbers, even for well documented drugs, so this always
#' attempts a PubChem-based CAS lookup as a fallback if `result$cas` is
#' still `NA` after the initial ChEMBL/PubChem resolution, using whichever
#' PubChem CID is available (or looking one up by SMILES if not).
#'
#' @param result A result list as returned by the ChEMBL/PubChem query
#'   functions.
#' @return The same list, with `cas` filled in if a PubChem lookup
#'   succeeded.
#' @noRd
fill_cas_via_pubchem <- function(result) {
  if (!is.na(result$cas)) {
    return(result)
  }

  cid <- result$pubchem_cid
  if (is.null(cid) || is.na(cid)) {
    cid <- get_pubchem_cid_by_smiles(result$smiles)
  }

  if (!is.na(cid)) {
    result$cas <- fetch_cas_by_pubchem_cid(cid)
  }

  result
}

#' Detect whether a string looks like a CAS registry number
#'
#' Internal helper. CAS numbers follow the pattern digits-digits-digit,
#' e.g. "58-08-2".
#'
#' @param x A character string.
#' @return `TRUE` if `x` matches the CAS number pattern, `FALSE` otherwise.
#' @noRd
is_cas_number <- function(x) {
  grepl("^[0-9]{2,7}-[0-9]{2}-[0-9]$", x)
}
