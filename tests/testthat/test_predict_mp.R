
# Inputs are illustrative round numbers, not literature values.

make_row <- function(MW, LogP, PSA, HBD, HBA, pKa_acid, pKa_base, Fu, LogD74,
                     Type, BDDCS_class) {
  data.frame(Drug = "test", MW = MW, LogP = LogP, PSA = PSA, HBD = HBD,
             HBA = HBA, pKa_acid = pKa_acid, pKa_base = pKa_base, Fu = Fu,
             LogD74 = LogD74, Type = Type, BDDCS_class = BDDCS_class,
             stringsAsFactors = FALSE)
}

acid <- make_row(180.16, 1.31, 63.6, 1, 3, 3.5, NA, 0.20,
                 1.31 + log10(1 / (1 + 10^(7.4 - 3.5))), "acid", 3)
base <- make_row(285.4, 2.9, 30, 1, 3, NA, 9.0, 0.60,
                 2.9 + log10(1 / (1 + 10^(9.0 - 7.4))), "base", 1)
neutral <- make_row(194.19, -0.07, 58.4, 0, 6, NA, NA, 0.70, -0.07, "neutral", 1)
amph1 <- make_row(350, 2, 80, 2, 5, 8.0, 4.0, 0.40, 1.9, "ampholyte", 2)
amph2 <- make_row(350, 2, 80, 2, 5, 5.0, 9.5, 0.40, 0.5, "ampholyte", 2)

test_that("mechanistic models match the reference implementation", {
  a <- predict_mp(acid)
  expect_equal(a$MuPu, 0.39818, tolerance = 1e-3)
  expect_equal(a$Phase_Distribution, 0.08207, tolerance = 1e-3)
  expect_equal(a$Log_PD, 0.03061, tolerance = 1e-3)
  expect_equal(a$Koshimichi, 0.0322, tolerance = 1e-3)
  expect_equal(a$Meskin_Lien, 0.44795, tolerance = 1e-3)

  b <- predict_mp(base)
  expect_equal(b$MuPu, 2.47484, tolerance = 1e-3)
  expect_equal(b$Phase_Distribution, 1.61263, tolerance = 1e-3)
  expect_equal(b$Koshimichi, 1.88466, tolerance = 1e-3)
  expect_equal(b$Meskin_Lien, 1.06194, tolerance = 1e-3)
})

test_that("neutral drugs have MuPu 1 and no Meskin & Lien value", {
  n <- predict_mp(neutral)
  expect_equal(n$MuPu, 1)
  expect_equal(n$Phase_Distribution, 0.70203, tolerance = 1e-3)
  expect_true(is.na(n$Meskin_Lien))
  expect_match(n$notes, "Meskin_Lien has no equation")
})

test_that("Meskin & Lien ampholyte rules (type 1 and type 2)", {
  expect_equal(predict_mp(amph1)$Meskin_Lien, 0.04648, tolerance = 1e-3)
  expect_equal(predict_mp(amph2)$Meskin_Lien, 0.44854, tolerance = 1e-3)
})

test_that("ML model returns the expected values", {
  skip_if_not_installed("xgboost")
  # raw model output -0.26878 / -0.04124; scale set by ML_TARGET_SCALE
  expect_equal(predict_mp(acid)$ML, back_transform_ml(-0.26878), tolerance = 1e-3)
  expect_equal(predict_mp(neutral)$ML, back_transform_ml(-0.04124), tolerance = 1e-3)
})

test_that("missing Fu leaves Fu-dependent methods NA but not Meskin & Lien", {
  no_fu <- acid
  no_fu$Fu <- NA
  r <- predict_mp(no_fu)
  expect_true(is.na(r$Phase_Distribution))
  expect_true(is.na(r$Log_PD))
  expect_true(is.na(r$Koshimichi))
  expect_false(is.na(r$Meskin_Lien))
  expect_match(r$notes, "Fu missing")
})

test_that("unknown ion class is not treated as neutral", {
  unk <- neutral
  unk$Type <- NA
  r <- predict_mp(unk)
  expect_true(is.na(r$MuPu))
  expect_true(is.na(r$Phase_Distribution))
})

test_that("bad input is rejected clearly", {
  expect_error(predict_mp(data.frame()), "at least one row")
  expect_error(predict_mp(acid[, c("MW", "LogP")]), "missing column")
})

test_that("Type labels and BDDCS formats are normalised", {
  expect_equal(normalise_ion_class("Amphoteric", NA, NA), "ampholyte")
  expect_equal(normalise_ion_class(NA, 4, NA), "acid")
  expect_true(is.na(normalise_ion_class(NA, NA, NA)))
  expect_equal(parse_bddcs(c("2", "Class 3", 4, NA, "n/a")), c(2L, 3L, 4L, NA, NA))
})
