#' Query ChEMBL by drug name
#'
#' Internal helper. Searches ChEMBL's molecule search endpoint for a name
#' match and extracts the parent structure's canonical SMILES.
#'
#' @param name A drug name.
#' @return A list matching the structure documented in [resolve_compound()],
#'   or `NULL` if no match was found.
#' @noRd
query_chembl_by_name <- function(name) {
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q=",
    utils::URLencode(name, reserved = TRUE)
  )
  parse_chembl_response(safe_get_json(url), name)
}

#' Query ChEMBL by CAS number
#'
#' Internal helper. ChEMBL does not index CAS numbers directly for full-text
#' search in the same way as names, so this searches on the CAS string as a
#' synonym/xref query.
#'
#' @param cas A CAS registry number.
#' @return A list matching the structure documented in [resolve_compound()],
#'   or `NULL` if no match was found.
#' @noRd
query_chembl_by_cas <- function(cas) {
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q=",
    utils::URLencode(cas, reserved = TRUE)
  )
  parse_chembl_response(safe_get_json(url), cas)
}

#' Parse a ChEMBL molecule search response
#'
#' Internal helper shared by the name- and CAS-based ChEMBL queries.
#'
#' @param resp Parsed JSON response (list), or `NULL` if the request failed.
#' @param query_input The original input string, used for the error/log
#'   message trail only.
#' @return A list matching [resolve_compound()]'s structure, or `NULL`.
#' @noRd
parse_chembl_response <- function(resp, query_input) {

  if (is.null(resp) || is.null(resp$molecules) || length(resp$molecules) == 0) {
    return(NULL)
  }

  mol <- select_best_chembl_match(resp$molecules, query_input)

  smiles <- tryCatch(mol$molecule_structures$canonical_smiles, error = function(e) NULL)
  if (is.null(smiles) || !nzchar(smiles)) {
    return(NULL)  # a "match" with no usable structure is not useful downstream
  }

  matched_name <- tryCatch(mol$pref_name, error = function(e) NULL)
  if (is.null(matched_name) || is.na(matched_name)) matched_name <- query_input

  chembl_id <- mol$molecule_chembl_id

  cas <- extract_cas_from_chembl_synonyms(mol)

  if (is.na(cas)) {
    # The molecule/search.json endpoint does not reliably populate
    # molecule_synonyms. Fall back to a direct fetch of the full molecule
    # record by ID, which does include synonyms per ChEMBL's schema.
    cas <- fetch_cas_by_chembl_id(chembl_id)
  }

  list(
    smiles       = smiles,
    cas          = cas,
    chembl_id    = chembl_id,
    pubchem_cid  = NA_character_,
    source       = "chembl",
    matched_name = matched_name
  )
}

#' Select the best-matching molecule from a ChEMBL search result list
#'
#' Internal helper used by [parse_chembl_response()]. ChEMBL's search
#' endpoint ranks results by text-search relevance, not exact match. To avoid silently resolving to the wrong
#' compound, an exact case insensitive match on `pref_name` is always
#' preferred over the top-ranked one, when one exists in the result
#' set.
#'
#' @param molecules The `molecules` list from a ChEMBL search response.
#' @param query_input The original search string.
#' @return The selected molecule record (a single list element).
#' @noRd
select_best_chembl_match <- function(molecules, query_input) {
  query_lower <- tolower(trimws(query_input))

  pref_names <- vapply(
    molecules,
    function(m) {
      pn <- tryCatch(m$pref_name, error = function(e) NULL)
      if (is.null(pn) || is.na(pn)) "" else tolower(pn)
    },
    character(1)
  )

  exact_idx <- which(pref_names == query_lower)

  if (length(exact_idx) > 0) {
    return(molecules[[exact_idx[1]]])
  }

  molecules[[1]]  # no exact match found; fall back to the top ranked
}

#' Fetch CAS number via a direct molecule by ID lookup
#'
#' Internal helper used as a fallback when [extract_cas_from_chembl_synonyms()]
#' finds nothing in a search-result molecule object (the search endpoint
#' does not reliably include molecule_synonyms). Queries the detailed
#' single-molecule endpoint instead, which does.
#'
#' @param chembl_id A ChEMBL molecule ID.
#' @return A CAS number string, or `NA` if none is listed.
#' @noRd
fetch_cas_by_chembl_id <- function(chembl_id) {
  if (is.null(chembl_id) || is.na(chembl_id)) return(NA_character_)

  url <- paste0("https://www.ebi.ac.uk/chembl/api/data/molecule/", chembl_id, ".json")
  resp <- safe_get_json(url)

  if (is.null(resp)) return(NA_character_)

  extract_cas_from_chembl_synonyms(resp)
}

#' Extract a CAS number from a ChEMBL molecule's synonym list, if present
#'
#' @param mol A single molecule record from a ChEMBL API response.
#' @return A CAS number string, or `NA` if none is listed.
#' @noRd
extract_cas_from_chembl_synonyms <- function(mol) {
  syns <- tryCatch(mol$molecule_synonyms, error = function(e) NULL)
  if (is.null(syns) || length(syns) == 0) return(NA_character_)

  cas_syn <- Filter(function(s) {
    identical(s$syn_type, "CAS_NUMBER") ||
      (!is.null(s$synonyms) && is_cas_number(s$synonyms))
  }, syns)

  if (length(cas_syn) == 0) return(NA_character_)
  cas_syn[[1]]$synonyms
}


# PHYSICAL PROPERTIES (used by extract_physical_properties())


#' Query ChEMBL for a compound's physical properties by ChEMBL ID
#'
#' Internal helper used by [extract_physical_properties()]. Fetches the
#' molecule record directly by ID and extracts the molecule_properties fields
#'
#' @param chembl_id A ChEMBL molecule ID.
#' @return A list with mw, logp, psa, hbd, hba (each numeric or NA).
#' @noRd
query_chembl_properties <- function(chembl_id) {
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/molecule/", chembl_id, ".json"
  )

  resp <- safe_get_json(url)

  if (is.null(resp) || is.null(resp$molecule_properties)) {
    return(empty_properties())
  }

  props <- resp$molecule_properties

  list(
    mw   = safe_numeric(props$full_mwt),
    logp = safe_numeric(props$alogp),
    psa  = safe_numeric(props$psa),
    hbd  = safe_integer(props$hbd),
    hba  = safe_integer(props$hba)
  )
}

#' Safely coerce a possibly NULL value to numeric
#' @noRd
safe_numeric <- function(x) {
  if (is.null(x)) return(NA_real_)
  v <- suppressWarnings(as.numeric(x))
  if (length(v) == 0) NA_real_ else v
}

#' Safely coerce a possibly NULL value to integer
#' @noRd
safe_integer <- function(x) {
  if (is.null(x)) return(NA_integer_)
  v <- suppressWarnings(as.integer(x))
  if (length(v) == 0) NA_integer_ else v
}
