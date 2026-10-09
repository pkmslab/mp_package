#' Look up a drug's BDDCS classification
#'
#' Matches a drug name against the Biopharmaceutics Drug Disposition
#' Classification System (BDDCS) reference table (Benet, Broccatelli &
#' Oprea, "BDDCS Applied to Over 900 Drugs", AAPS J. 2011;13(4):519-547,
#' \url{https://pmc.ncbi.nlm.nih.gov/articles/PMC3231854/}), covering 927
#' drugs.
#'
#' @param name A drug name to look up. Matching is case insensitive and
#'   exact  against the reference table's generic names. Some
#'   reference entries list two synonyms joined by `"; "` (e.g.
#'   `"Acetaminophen; Paracetamol"`)  these are matched if `name` equals
#'   either half.
#'
#' @return An integer BDDCS class (1-4), or:
#' \itemize{
#'   \item `0` if the drug is present in the reference table but the
#'     original source left it unclassified (11 drugs fall in this
#'     category, e.g. Amphetamine, Vitamin C)
#'   \item `NA` if the name has no match in the reference table at all ,
#'     this includes drugs known only by an alternate name not listed as
#'     a synonym, which should be
#'     entered manually
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' lookup_bddcs("Aspirin")        # matches via the "Acetylsalicylic
#'                                 # Acid; Aspirin" synonym entry
#' lookup_bddcs("Methyldopa")     # exact match, class 3
#' }
lookup_bddcs <- function(name) {

  if (missing(name) || !is.character(name) || length(name) != 1 || !nzchar(trimws(name))) {
    stop("`name` must be a single non-empty character string.", call. = FALSE)
  }

  table <- load_bddcs_table()
  name_lower <- tolower(trimws(name))

  match_idx <- which(table$name_lower == name_lower)

  if (length(match_idx) == 0) {
    return(NA_integer_)
  }

  table$bddcs_class[match_idx[1]]
}

#' Load and cache the BDDCS reference table
#'
#' Internal helper. Reads the bundled CSV once per session, splits any
#' semicolon-joined synonym entries into separate matchable rows, and
#' caches the result so repeated calls to [lookup_bddcs()] don't re-read
#' the file from disk.
#'
#' @return A data.frame with columns `name_lower` and `bddcs_class`.
#' @noRd
load_bddcs_table <- function() {
  if (!is.null(.bddcs_cache$table)) {
    return(.bddcs_cache$table)
  }

  csv_path <- system.file("extdata", "bddcs_lookup.csv", package = "mp.predict")

  if (!nzchar(csv_path)) {
    stop(
      "Could not locate bddcs_lookup.csv inside the installed package. ",
      "Check that inst/extdata/ is present.",
      call. = FALSE
    )
  }

  raw <- utils::read.csv(csv_path, stringsAsFactors = FALSE)

  expanded <- do.call(rbind, lapply(seq_len(nrow(raw)), function(i) {
    synonyms <- trimws(strsplit(raw$generic_name[i], ";")[[1]])
    data.frame(
      name_lower  = tolower(synonyms),
      bddcs_class = raw$bddcs_class[i],
      stringsAsFactors = FALSE
    )
  }))

  .bddcs_cache$table <- expanded
  expanded
}

# Session-level cache environment for the BDDCS table
.bddcs_cache <- new.env(parent = emptyenv())
