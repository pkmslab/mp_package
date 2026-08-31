#' Call pkasolver on a single SMILES string
#'
#' Internal helper used by [extract_pka()]. Activates the target conda
#' environment, imports pkasolver and RDKit, and runs the microstate pKa
#' prediction for one compound.
#'
#' @param smiles A SMILES string.
#' @param envname Conda environment name.
#' @return A list with `pka_values` (numeric vector) and `n_sites` (integer).
#' @noRd
call_pkasolver <- function(smiles, envname = "pka_prediction") {

  reticulate::use_condaenv(envname, required = TRUE)

  rdkit_chem      <- reticulate::import("rdkit.Chem", delay_load = FALSE)
  pkasolver_query <- reticulate::import("pkasolver.query", delay_load = FALSE)

  mol <- rdkit_chem$MolFromSmiles(smiles)

  if (is.null(mol)) {
    warning("RDKit could not parse SMILES: '", smiles, "'. Returning no pKa values.", call. = FALSE)
    return(list(pka_values = numeric(0), n_sites = 0L))
  }

  states <- tryCatch(
    pkasolver_query$calculate_microstate_pka_values(mol, only_dimorphite = FALSE),
    error = function(e) {
      warning("pkasolver failed on '", smiles, "': ", conditionMessage(e), call. = FALSE)
      list()
    }
  )

  if (length(states) == 0) {
    return(list(pka_values = numeric(0), n_sites = 0L))
  }

  pka_values <- vapply(states, function(s) round(as.numeric(s$pka), 3), numeric(1))

  list(pka_values = pka_values, n_sites = length(pka_values))
}
