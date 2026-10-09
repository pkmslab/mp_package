#' Safely GET a URL and parse it as JSON
#'
#' Internal helper shared across all API query functions. Wraps the request
#' in error handling so that network failures, timeouts, or non-200
#' responses return `NULL` instead of crashing the calling function
#' letting [resolve_compound()] cleanly fall through to the next data
#' source.
#'
#' @param url The full request URL.
#' @param timeout_sec Request timeout in seconds. Defaults to 10.
#' @return The parsed JSON response as a list, or `NULL` on any failure
#'   (network error, timeout, non-200 status, or unparseable body).
#' @noRd
safe_get_json <- function(url, timeout_sec = 10) {
  tryCatch({
    resp <- httr::GET(url, httr::timeout(timeout_sec))

    if (httr::status_code(resp) != 200) {
      return(NULL)
    }

    body_text <- httr::content(resp, as = "text", encoding = "UTF-8")
    if (!nzchar(body_text)) return(NULL)

    jsonlite::fromJSON(body_text, simplifyVector = FALSE)
  }, error = function(e) {
    NULL
  })
}
