#' Predict milk to plasma (M/P) ratios
#'
#' Takes the feature table from [extract_mp_features()] (after any manual
#' edits) and returns one M/P prediction per method:
#'
#' * `ML` - XGBoost model trained on the M/P dataset
#' * `Phase_Distribution` - phase distribution model
#' * `Log_PD` - log-transformed phase distribution model
#' * `Koshimichi` - clearance-ratio model (Cl_sec / Cl_re)
#' * `Meskin_Lien` - Meskin & Lien regression (acids / bases)
#'
#' A method that cannot be computed from the available inputs returns `NA`
#' and the reason is written to the `notes` column. Nothing is guessed.
#'
#' @param table A data frame from [extract_mp_features()], one row per drug.
#'   Needs the columns `MW`, `LogP`, `PSA`, `HBD`, `HBA`, `pKa_acid`,
#'   `pKa_base`, `Fu`, `Type`, `LogD74` and `BDDCS_class`. Missing values
#'   (`NA`) are allowed.
#' @return A data frame with `Drug`, `MuPu`, the five predictions and `notes`.
#' @examples
#' \dontrun{
#' feats <- extract_mp_features("aspirin")
#' predict_mp(feats)
#' }
#' @export
predict_mp <- function(table) {
  if (!is.data.frame(table) || nrow(table) == 0) {
    stop("`table` must be a data frame with at least one row, ",
         "as returned by extract_mp_features().", call. = FALSE)
  }
  required <- c("MW", "LogP", "PSA", "HBD", "HBA", "pKa_acid", "pKa_base",
                "Fu", "Type", "LogD74", "BDDCS_class")
  missing_cols <- setdiff(required, names(table))
  if (length(missing_cols) > 0) {
    stop("`table` is missing column(s): ", paste(missing_cols, collapse = ", "),
         ". Use the output of extract_mp_features().", call. = FALSE)
  }

  num <- function(x) suppressWarnings(as.numeric(x))
  x <- data.frame(
    MW = num(table$MW), LogP = num(table$LogP), PSA = num(table$PSA),
    HBD = num(table$HBD), HBA = num(table$HBA),
    pKa_acid = num(table$pKa_acid), pKa_base = num(table$pKa_base),
    Fu = num(table$Fu), LogD74 = num(table$LogD74),
    stringsAsFactors = FALSE
  )
  x$Type <- normalise_ion_class(table$Type, x$pKa_acid, x$pKa_base)
  x$BDDCS <- parse_bddcs(table$BDDCS_class)

  # Fu outside (0, 1] is a data error (and Fu = 0 breaks the log models)
  bad_fu <- !is.na(x$Fu) & (x$Fu <= 0 | x$Fu > 1)
  x$Fu[bad_fu] <- NA_real_

  mech <- mechanistic_variables(x)
  x$MuPu <- mech$MuPu

  out <- data.frame(
    Drug = if ("Drug" %in% names(table)) as.character(table$Drug) else NA_character_,
    MuPu = mech$MuPu,
    ML = predict_ml(x),
    Phase_Distribution = mp_phase_distribution(x, mech),
    Log_PD = mp_log_pd(x, mech),
    Koshimichi = mp_koshimichi(x, mech),
    Meskin_Lien = mp_meskin_lien(x),
    stringsAsFactors = FALSE
  )
  out$notes <- build_notes(x, out, bad_fu)
  rownames(out) <- NULL
  out
}

# CONSTANTS

# Scale the XGBoost model was trained on: log10(M/P)
#   "log10" -> prediction = 10^output      "ln" -> prediction = exp(output)
ML_TARGET_SCALE <- "log10"

PH_PLASMA <- 7.4
PH_MILK   <- 7.0

# INPUT CLEANING

#' Standardize the ion class to neutral / acid / base / ampholyte.
#' Falls back on which pKa values are present when `Type` is missing or not
#' recognized. A drug with no pKa and no usable Type is left NA (unknown),
#' not "neutral".
#' @noRd
normalise_ion_class <- function(type, pka_acid, pka_base) {
  t <- tolower(trimws(as.character(type)))
  t[t %in% c("amphoteric", "zwitterion", "zwitterionic", "ampholytic")] <- "ampholyte"
  t[t %in% c("acidic")] <- "acid"
  t[t %in% c("basic")] <- "base"
  t[!(t %in% c("neutral", "acid", "base", "ampholyte"))] <- NA_character_

  from_pka <- ifelse(!is.na(pka_acid) & !is.na(pka_base), "ampholyte",
                     ifelse(!is.na(pka_acid), "acid",
                            ifelse(!is.na(pka_base), "base", NA_character_)))
  ifelse(is.na(t), from_pka, t)
}

#' "1", "Class 2", 3, "BDDCS 4" -> 1..4; anything else -> NA
#' @noRd
parse_bddcs <- function(x) {
  d <- suppressWarnings(as.integer(sub(".*?([1-4]).*", "\\1", as.character(x))))
  d[is.na(x) | !grepl("[1-4]", as.character(x))] <- NA_integer_
  d
}


# SHARED MECHANISTIC VARIABLES


#' Fraction un-ionized at a given pH. One formula for every ion class:
#' f = 1 / (1 + 10^(pH - pKa_acid) + 10^(pKa_base - pH)); a missing pKa
#' contributes nothing, so it reduces to the acid only or base only formula.
#' Returns NA when the drug's ion class is unknown (no pKa information).
#' @noRd
fraction_unionised <- function(pka_acid, pka_base, ion_class, pH) {
  acid_term <- ifelse(is.na(pka_acid), 0, 10^(pH - pka_acid))
  base_term <- ifelse(is.na(pka_base), 0, 10^(pka_base - pH))
  f <- 1 / (1 + acid_term + base_term)
  f[is.na(ion_class)] <- NA_real_
  f
}

#' @noRd
mechanistic_variables <- function(x) {
  fu_p <- x$Fu
  fun_p <- fraction_unionised(x$pKa_acid, x$pKa_base, x$Type, PH_PLASMA)
  fun_m <- fraction_unionised(x$pKa_acid, x$pKa_base, x$Type, PH_MILK)

  fu_m <- fu_p^0.448 / (0.000694^0.448 + fu_p^0.448)
  mupu <- fun_p / fun_m

  d70      <- 10^x$LogD74 * fun_m / fun_p
  logd70   <- -0.88 + 1.29 * log10(d70)
  d_milk   <- 10^logd70

  log_cl_sec <- -3.912 - 0.015 * x$PSA + 3.367 * log10(x$MW) -
    0.164 * log10(10^x$LogP / 10^x$LogD74)
  log_cl_re  <- 2.793 + 0.179 * x$LogP - 0.132 * x$HBD
  cl_ratio   <- 10^log_cl_sec / 10^log_cl_re

  f_m <- 1 / (0.955 / fu_m + 0.045 * d_milk)

  list(Fu_m = fu_m, fun_p = fun_p, fun_m = fun_m, MuPu = mupu,
       D_milk = d_milk, Cl_ratio = cl_ratio, f_m = f_m)
}


# THE FOUR MECHANISTIC MODELS


#' MP = Fu_p x Mu/Pu x (0.955/Fu_m + 0.045 x D_milk)
#' @noRd
mp_phase_distribution <- function(x, mech) {
  x$Fu * mech$MuPu * (0.955 / mech$Fu_m + 0.045 * mech$D_milk)
}

#' ln(MP) = 0.03 + 2.28 ln(Mu/Pu) + 0.89 ln(Fu_p) + 0.51 ln(0.955/Fu_m + 0.045 D_milk)
#' @noRd
mp_log_pd <- function(x, mech) {
  ln_mp <- 0.03 + 2.28 * log(mech$MuPu) + 0.89 * log(x$Fu) +
    0.51 * log(0.955 / mech$Fu_m + 0.045 * mech$D_milk)
  exp(ln_mp)
}

#' MP = (Cl_sec / Cl_re) x (Fu_p / f_m)
#' @noRd
mp_koshimichi <- function(x, mech) {
  mech$Cl_ratio * (x$Fu / mech$f_m)
}

#' Meskin & Lien. Acids: 10^(2.068 - 0.162 sqrt(MW) - 0.185 LogP).
#' Bases: 10^(0.265 - 0.153 LogP - 0.128 (7.4 - pKa)).
#' Ampholytes are assigned to the acid or base equation by whichever pKa is
#' closer to 7 (the pH of breast milk):
#'   type 1 (acidic pKa > basic pKa): closer to 7 decides, basic pKa used
#'   type 2 (basic pKa > acidic pKa): basic closer to 7 -> acid equation;
#'                                    acidic closer to 7 -> base equation
#'                                    using the ACIDIC pKa
#' Neutral drugs (no ionisable group) have no Meskin & Lien equation: NA.
#' @noRd
mp_meskin_lien <- function(x) {
  n <- nrow(x)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    cls <- x$Type[i]
    if (is.na(cls) || cls == "neutral") next

    eq <- NA_character_
    pka <- NA_real_
    pa <- x$pKa_acid[i]
    pb <- x$pKa_base[i]

    if (cls == "acid") {
      eq <- "acid"
    } else if (cls == "base") {
      eq <- "base"; pka <- pb
    } else if (cls == "ampholyte") {
      if (is.na(pa) && is.na(pb)) {
        next
      } else if (is.na(pb)) {            # only an acidic site left
        eq <- "acid"
      } else if (is.na(pa)) {            # only a basic site left
        eq <- "base"; pka <- pb
      } else {
        acid_closer <- abs(pa - 7) < abs(pb - 7)
        if (pa > pb) {                   # type 1
          if (acid_closer) { eq <- "acid" } else { eq <- "base"; pka <- pb }
        } else {                         # type 2
          if (acid_closer) { eq <- "base"; pka <- pa } else { eq <- "acid" }
        }
      }
    }

    if (eq == "acid") {
      out[i] <- 10^(2.068 - 0.162 * sqrt(x$MW[i]) - 0.185 * x$LogP[i])
    } else if (eq == "base") {
      out[i] <- 10^(0.265 - 0.153 * x$LogP[i] - 0.128 * (7.4 - pka))
    }
  }
  out
}

# MACHINE LEARNING MODEL


ML_FEATURES <- c("MW", "PSA", "HBD", "HBA", "LogP", "Fu", "pKa_acid", "pKa_base",
                 "MuPu", "LogD74", "BDDCS_1", "BDDCS_2", "BDDCS_3", "BDDCS_4",
                 "ion_class_neutral", "ion_class_acid", "ion_class_base",
                 "ion_class_ampholyte")

#' Build the 18-column design matrix in the exact order the model was
#' trained on. Unknown BDDCS / ion class is left NA (XGBoost sends it down
#' its default branch) rather than filled with zeros, which would claim a
#' class that was never established.
#' @noRd
ml_design_matrix <- function(x) {
  m <- cbind(
    MW = x$MW, PSA = x$PSA, HBD = x$HBD, HBA = x$HBA, LogP = x$LogP,
    Fu = x$Fu, pKa_acid = x$pKa_acid, pKa_base = x$pKa_base,
    MuPu = x$MuPu, LogD74 = x$LogD74,
    BDDCS_1 = as.numeric(x$BDDCS == 1), BDDCS_2 = as.numeric(x$BDDCS == 2),
    BDDCS_3 = as.numeric(x$BDDCS == 3), BDDCS_4 = as.numeric(x$BDDCS == 4),
    ion_class_neutral   = as.numeric(x$Type == "neutral"),
    ion_class_acid      = as.numeric(x$Type == "acid"),
    ion_class_base      = as.numeric(x$Type == "base"),
    ion_class_ampholyte = as.numeric(x$Type == "ampholyte")
  )
  # `==` on NA gives NA, which is what we want for unknown classes
  m[, ML_FEATURES, drop = FALSE]
}

#' @noRd
load_ml_model <- function() {
  path <- system.file("extdata", "final_MP_model.json", package = "mp.predict")
  if (!nzchar(path)) {
    stop("Model file not found in the installed package ",
         "(inst/extdata/final_MP_model.json).", call. = FALSE)
  }
  xgboost::xgb.load(path)
}

#' @noRd
predict_ml <- function(x) {
  model <- load_ml_model()
  raw <- as.numeric(stats::predict(model, ml_design_matrix(x)))
  back_transform_ml(raw)
}

#' @noRd
back_transform_ml <- function(raw) {
  switch(ML_TARGET_SCALE,
         log10 = 10^raw,
         ln    = exp(raw),
         stop("ML_TARGET_SCALE must be 'log10' or 'ln'.", call. = FALSE))
}

# Notes: say why a method came back NA

#' @noRd
build_notes <- function(x, out, bad_fu) {
  vapply(seq_len(nrow(x)), function(i) {
    n <- character(0)
    if (bad_fu[i]) n <- c(n, "Fu outside (0, 1] ignored")
    if (is.na(x$Fu[i])) n <- c(n, "Fu missing: Phase_Distribution, Log_PD and Koshimichi not computed")
    if (is.na(x$Type[i])) n <- c(n, "no pKa or ion class: mechanistic models not computed")
    if (!is.na(x$Type[i]) && x$Type[i] == "neutral") n <- c(n, "neutral drug: Meskin_Lien has no equation")
    if (any(is.na(x[i, c("MW", "LogP", "PSA", "HBD", "LogD74")]))) {
      n <- c(n, "MW/LogP/PSA/HBD/LogD74 missing: some methods not computed")
    }
    if (is.na(x$BDDCS[i])) n <- c(n, "BDDCS class unknown (ML uses its default branch)")
    if (length(n) == 0) NA_character_ else paste(n, collapse = "; ")
  }, character(1))
}
