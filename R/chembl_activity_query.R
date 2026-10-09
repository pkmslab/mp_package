#' Query ChEMBL activity data for direct Fu measurements
#'
#' Internal helper used by [extract_fu()]. Ports the exact filtering logic
#' used to build the training dataset: human-only, plasma-only assays,
#' animal keywords excluded, values restricted to the 0-1 range.
#'
#' @param chembl_id A ChEMBL molecule ID.
#' @return A list with `values` (numeric vector of valid Fu values),
#'   `n_values` (integer), and `flagged` (TRUE if any raw value was exactly
#'   0 or 1).
#' @noRd
query_chembl_fu <- function(chembl_id) {
  activities <- fetch_chembl_activities(chembl_id, standard_type = "Fu")

  if (is.null(activities) || nrow(activities) == 0) {
    return(list(values = numeric(0), n_values = 0L, flagged = FALSE))
  }

  filtered <- filter_human_plasma_activities(activities, exclude_tissue = TRUE)
  filtered <- filtered[filtered$standard_value >= 0 & filtered$standard_value <= 1, , drop = FALSE]

  if (nrow(filtered) == 0) {
    return(list(values = numeric(0), n_values = 0L, flagged = FALSE))
  }

  list(
    values   = filtered$standard_value,
    n_values = nrow(filtered),
    flagged  = any(filtered$standard_value == 0 | filtered$standard_value == 1)
  )
}

#' Query ChEMBL activity data for plasma protein binding (PPB), as a Fu fallback
#'
#' Internal helper used by [extract_fu()] when no direct Fu measurement is
#' available. Same human/plasma filtering as [query_chembl_fu()], but values
#' are restricted to the 0-100 range (percent bound) and converted to Fu via
#' `Fu = 1 - (PPB / 100)`.
#'
#' @param chembl_id A ChEMBL molecule ID.
#' @return A list with `fu_values` (numeric vector, already converted from
#'   PPB) and `n_values` (integer).
#' @noRd
query_chembl_ppb <- function(chembl_id) {
  activities <- fetch_chembl_activities(chembl_id, standard_type = "PPB")

  if (is.null(activities) || nrow(activities) == 0) {
    return(list(fu_values = numeric(0), n_values = 0L))
  }

  filtered <- filter_human_plasma_activities(activities, exclude_tissue = FALSE)
  filtered <- filtered[filtered$standard_value >= 0 & filtered$standard_value <= 100, , drop = FALSE]

  if (nrow(filtered) == 0) {
    return(list(fu_values = numeric(0), n_values = 0L))
  }

  list(
    fu_values = 1 - (filtered$standard_value / 100),
    n_values  = nrow(filtered)
  )
}

#' Fetch raw activity records from ChEMBL for a given standard_type
#'
#' Internal helper shared by [query_chembl_fu()] and [query_chembl_ppb()].
#'
#' @param chembl_id A ChEMBL molecule ID.
#' @param standard_type ChEMBL's `standard_type` filter (e.g. `"Fu"`,
#'   `"PPB"`).
#' @return A data.frame of activity records with at least
#'   `target_organism`, `assay_description`, `standard_value` columns, or
#'   `NULL` if the request failed or returned nothing.
#' @noRd
fetch_chembl_activities <- function(chembl_id, standard_type) {
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id=",
    chembl_id, "&standard_type=", standard_type, "&limit=100"
  )

  resp <- safe_get_json(url)
  if (is.null(resp) || is.null(resp$activities) || length(resp$activities) == 0) {
    return(NULL)
  }

  activities_to_dataframe(resp$activities)
}

#' Convert ChEMBL's list of lists activity records into a data.frame
#'
#' @param activities A list of activity records, as parsed from ChEMBL's
#'   JSON response.
#' @return A data.frame with columns target_organism, assay_description,
#'   standard_value (numeric).
#' @noRd
activities_to_dataframe <- function(activities) {
  data.frame(
    target_organism = vapply(
      activities, function(a) if (is.null(a$target_organism)) NA_character_ else a$target_organism,
      character(1)
    ),
    assay_description = vapply(
      activities, function(a) if (is.null(a$assay_description)) "" else a$assay_description,
      character(1)
    ),
    standard_value = vapply(
      activities,
      function(a) {
        v <- suppressWarnings(as.numeric(a$standard_value))
        if (length(v) == 0 || is.null(a$standard_value)) NA_real_ else v
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )
}

#' Apply the human-only, plasma-only, animal/tissue-excluded filter
#'
#' Internal helper shared by [query_chembl_fu()] and [query_chembl_ppb()].
#' Ports the exact filtering rules used to build the training dataset.
#'
#' @param activities A data.frame as returned by [activities_to_dataframe()].
#' @param exclude_tissue If `TRUE`, also exclude records whose
#'   assay_description mentions "tissue" (used for Fu, not PPB, matching
#'   the original dataset logic).
#' @return The filtered data.frame.
#' @noRd
filter_human_plasma_activities <- function(activities, exclude_tissue = TRUE) {

  desc_lower <- tolower(activities$assay_description)

  is_human <- is.na(activities$target_organism) | activities$target_organism == "Homo sapiens"

  animal_pattern <- "\\brat\\b|\\bmouse\\b|\\bmice\\b|\\banimal\\b|\\brabbit\\b|\\bdog\\b|\\bmonkey\\b"
  no_animal_keywords <- !grepl(animal_pattern, desc_lower)

  no_tissue <- if (exclude_tissue) !grepl("tissue", desc_lower) else TRUE

  mentions_plasma_or_human <- grepl("plasma|human", desc_lower)

  has_value <- !is.na(activities$standard_value)

  keep <- is_human & no_animal_keywords & no_tissue & mentions_plasma_or_human & has_value

  activities[keep, , drop = FALSE]
}
