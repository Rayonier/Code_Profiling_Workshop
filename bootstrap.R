# One-time setup for the profiling workshop.
# Run from the workshop folder:  source("bootstrap.R")
# Installs any missing packages from CRAN, then checks that C++ compiles.

need <- c("Rcpp", "profvis", "bench", "rmarkdown", "knitr", "ggplot2", "pkgbuild")

miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  message("Installing: ", paste(miss, collapse = ", "))
  install.packages(miss, repos = "https://cloud.r-project.org")
  still <- miss[!vapply(miss, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still)) stop("These packages failed to install: ", paste(still, collapse = ", "))
}

if (!pkgbuild::has_build_tools(debug = FALSE)) {
  stop("No C++ toolchain found. On Windows, install the Rtools version that matches ",
       "your R (", R.version$major, ".", R.version$minor, ") from ",
       "https://cran.r-project.org/bin/windows/Rtools/ and restart R.")
}

message("Compiling a test C++ function (first time can take a minute)...")
Rcpp::cppFunction("int add_one(int x) { return x + 1; }")
if (add_one(1L) != 2L) stop("The test C++ function compiled but returned the wrong answer.")

if (!file.exists(file.path("data", "trees.csv"))) {
  stop("data/trees.csv not found. Set the working directory to the workshop folder ",
       "(RStudio: Session > Set Working Directory > To Source File Location) and re-run.")
}

cat("Setup OK -- you're ready.\n")
