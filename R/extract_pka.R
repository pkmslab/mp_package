#' Extract acidic and basic pKa values for a compound
#'
#' Predicts ionizable sites for a molecule using the pkasolver Python
#' package (Mayr et al.), then classifies each site as acidic or basic
#' (via charge comparison between protonated/deprotonated structures) and
#' identifies its functional group (via SMARTS substructure matching).
#' Only sites with pKa in the physiologically relevant range 3-11 AND a
#' recognized, approved functional group are kept. If multiple valid sites
#' remain, the lowest pKa among acidic sites and the highest pKa among
#' basic sites are returned as the final values.
#'
#' @param smiles A SMILES string for the (parent, salt-stripped) compound.
#' @param envname Name of the conda environment created by
#'   [setup_pkasolver()]. Defaults to `"pka_prediction"`.
#'
#' @return A list with components:
#' \describe{
#'   \item{pka_acid}{The lowest pKa among valid acidic sites, or `NA` if
#'     none found.}
#'   \item{pka_base}{The highest pKa among valid basic sites, or `NA` if
#'     none found.}
#'   \item{fg_acid}{Name of the functional group at the chosen acidic site
#'     (e.g. `"carboxylate"`), or `NA`.}
#'   \item{fg_base}{Name of the functional group at the chosen basic site,
#'     or `NA`.}
#'   \item{all_sites}{A data.frame of every ionizable site pkasolver found,
#'     including ones excluded by the range or functional-group filters --
#'     useful for double-checking a specific drug's classification.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' setup_pkasolver()  # once per machine
#' extract_pka("CN1C=NC2=C1C(=O)N(C)C(=O)N2C")  # caffeine
#' }
extract_pka <- function(smiles, envname = "pka_prediction") {

  if (missing(smiles) || !is.character(smiles) || length(smiles) != 1 || !nzchar(smiles)) {
    stop("`smiles` must be a single non-empty character string.", call. = FALSE)
  }

  if (!check_pkasolver_installed(envname)) {
    stop(
      "pkasolver is not set up in conda environment '", envname, "'.\n",
      "Run setup_pkasolver() once before calling extract_pka().",
      call. = FALSE
    )
  }

  result <- call_pkasolver(smiles, envname = envname)

  if (is.na(result$pka_acid) && is.na(result$pka_base)) {
    message("No valid ionizable sites found (after range/functional-group filtering); both pKa_acid and pKa_base returned as NA.")
  }

  result
}
