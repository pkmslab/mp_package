#' Fetch CAS number(s) for a PubChem CID
#'
#' Internal helper used by [resolve_compound()] as the primary CAS source.
#'
#' Uses PubChem's `/synonyms` endpoint (ranked/priority order) rather than
#' `/xrefs/RN` (all registry numbers ever assigned, unordered) , `xrefs/RN`
#' was found, via a real response for Aspirin's CID, to list a deprecated
#' CAS number ("11126-35-5") before the current one ("50-78-2"), with no
#' way to distinguish them. The synonyms endpoint's ranking reliably
#' surfaces the primary/current CAS number first in practice.
#'
#' @param cid A PubChem CID.
#' @return The first CAS-formatted entry in the compound's ranked synonym
#'   list, or `NA` if none found.
#' @noRd
fetch_cas_by_pubchem_cid <- function(cid) {
  if (is.null(cid) || is.na(cid)) return(NA_character_)

  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/synonyms/JSON"
  )
  resp <- safe_get_json(url)

  synonyms <- tryCatch(
    resp$InformationList$Information[[1]]$Synonym,
    error = function(e) NULL
  )

  if (is.null(synonyms) || length(synonyms) == 0) return(NA_character_)

  synonyms <- unlist(synonyms)
  cas_candidates <- Filter(is_cas_number, synonyms)
  if (length(cas_candidates) == 0) return(NA_character_)

  cas_candidates[[1]]
}

#' Look up a PubChem CID from a SMILES string
#'
#' Internal helper used by [resolve_compound()] to get a CID for CAS
#' lookup purposes, regardless of whether the compound was originally
#' resolved via ChEMBL or PubChem.
#'
#' @param smiles A SMILES string.
#' @return A PubChem CID (character), or `NA` if not found.
#' @noRd
get_pubchem_cid_by_smiles <- function(smiles) {
  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
    utils::URLencode(smiles, reserved = TRUE),
    "/cids/JSON"
  )
  resp <- safe_get_json(url)
  cid <- tryCatch(resp$IdentifierList$CID[[1]], error = function(e) NULL)
  if (is.null(cid)) return(NA_character_)
  as.character(cid)
}

#' Query PubChem by drug name
#'
#' Internal helper. Uses PubChem's PUG REST API to resolve a name to a CID
#' and canonical SMILES.
#'
#' @param name A drug name.
#' @return A list matching [resolve_compound()]'s structure, or `NULL` if no
#'   match was found.
#' @noRd
query_pubchem_by_name <- function(name) {
  cid_url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
    utils::URLencode(name, reserved = TRUE),
    "/cids/JSON"
  )
  cid_resp <- safe_get_json(cid_url)
  cid <- tryCatch(cid_resp$IdentifierList$CID[[1]], error = function(e) NULL)
  if (is.null(cid)) return(NULL)

  build_pubchem_result(cid, matched_name = name)
}

#' Query PubChem by CAS number
#'
#' Internal helper. Uses PubChem's xref endpoint to resolve a CAS registry
#' number to a CID.
#'
#' @param cas A CAS registry number.
#' @return A list matching [resolve_compound()]'s structure, or `NULL` if no
#'   match was found.
#' @noRd
query_pubchem_by_cas <- function(cas) {
  cid_url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/RN/",
    utils::URLencode(cas, reserved = TRUE),
    "/cids/JSON"
  )
  cid_resp <- safe_get_json(cid_url)
  cid <- tryCatch(cid_resp$IdentifierList$CID[[1]], error = function(e) NULL)
  if (is.null(cid)) return(NULL)

  build_pubchem_result(cid, matched_name = cas, known_cas = cas)
}

#' Fetch canonical SMILES for a PubChem CID and assemble the result list
#'
#' Internal helper shared by the name- and CAS-based PubChem queries.
#'
#' @param cid A PubChem compound ID.
#' @param matched_name The identifier used in the original query, for
#'   logging/reporting.
#' @param known_cas If the CAS number is already known (i.e. it was the
#'   query input), pass it through directly rather than re deriving it.
#' @return A list matching [resolve_compound()]'s structure, or `NULL` if the
#'   SMILES lookup failed.
#' @noRd
build_pubchem_result <- function(cid, matched_name, known_cas = NA_character_) {
  smiles_url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid,
    "/property/CanonicalSMILES/JSON"
  )
  smiles_resp <- safe_get_json(smiles_url)
  smiles <- tryCatch(
    smiles_resp$PropertyTable$Properties[[1]]$CanonicalSMILES,
    error = function(e) NULL
  )

  if (is.null(smiles) || !nzchar(smiles)) return(NULL)

  list(
    smiles       = smiles,
    cas          = known_cas,
    chembl_id    = NA_character_,
    pubchem_cid  = as.character(cid),
    source       = "pubchem",
    matched_name = matched_name
  )
}


# PHYSICAL PROPERTIES (used by extract_physical_properties())


#' Query PubChem for a compound's physical properties by SMILES
#'
#' Internal helper used by [extract_physical_properties()] as a fallback
#' when ChEMBL is missing some or all properties. Queries PubChem's PUG
#' REST API directly by structure (SMILES), so it works regardless of
#' whether a PubChem CID was already resolved.
#'
#' @param smiles A SMILES string.
#' @return A list with mw, logp, psa, hbd, hba (each numeric or NA).
#' @noRd
query_pubchem_properties_by_smiles <- function(smiles) {
  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/",
    utils::URLencode(smiles, reserved = TRUE),
    "/property/MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount/JSON"
  )

  resp <- safe_get_json(url)

  props <- tryCatch(resp$PropertyTable$Properties[[1]], error = function(e) NULL)

  if (is.null(props)) {
    return(empty_properties())
  }

  list(
    mw   = safe_numeric(props$MolecularWeight),
    logp = safe_numeric(props$XLogP),
    psa  = safe_numeric(props$TPSA),
    hbd  = safe_integer(props$HBondDonorCount),
    hba  = safe_integer(props$HBondAcceptorCount)
  )
}
