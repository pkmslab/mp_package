.onLoad <- function(libname, pkgname) {
  if (!nzchar(Sys.getenv("KMP_DUPLICATE_LIB_OK"))) {
    Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")
  }
}

# pkasolver (via reticulate -> torch) and xgboost each bundle their own
# OpenMP runtime. Loading both into one R process aborts the session on
# macOS ("R Session Aborted"). This setting lets the two coexist. It is set
# before any Python is started, and only if the user has not set it already.
