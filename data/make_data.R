# Generates the workshop's simulated loblolly pine plantation plots.
# Deterministic: rerunning gives identical files. Run from the workshop folder:
#   Rscript data/make_data.R
#
# What is made up vs. published:
#   * Made up (simulated stands): number of plots and trees, stand age, site
#     index, planting density, the DBH distribution, and crown ratios.
#   * Published: tree heights come from the FVS Southern variant Curtis-Arney
#     height-diameter equation for loblolly pine {4.1.1}, rescaled so each
#     plot's top height matches the published Carmean et al. (1989) site curve
#     {4.6.1.1} for its age and site index. Both are from:
#     FVS Staff (2008, revised 2026). Southern (SN) Variant Overview - Forest
#     Vegetation Simulator. USDA Forest Service.

set.seed(20261028)

N_PLOTS    <- 2000          # plots to simulate
TREES_MIN  <- 150           # tree records per plot
TREES_MAX  <- 250

# Published loblolly pine (LP) coefficients
CA <- c(P2 = 243.860648, P3 = 4.28460566, P4 = -0.47130185)          # Table 4.1.1
CM <- c(c1 = 1.1421, c2 = 1.0042, c3 = -0.0374, c4 = 0.7632, c5 = 0.0358)  # Table 4.6.1.1

curtis_arney <- function(dbh) 4.5 + CA[["P2"]] * exp(-CA[["P3"]] * dbh^CA[["P4"]])
carmean <- function(age, si) {
  CM[["c1"]] * si^CM[["c2"]] * (1 - exp(CM[["c3"]] * age))^(CM[["c4"]] * si^CM[["c5"]])
}

plot_ids <- sprintf("P%04d", seq_len(N_PLOTS))
plots <- data.frame(
  PlotID = plot_ids,
  Age    = sample(8:15, N_PLOTS, replace = TRUE),
  SI     = round(runif(N_PLOTS, 55, 85))          # site index, base age 50 (ft)
)

make_trees <- function(id, age, si) {
  n       <- sample(TREES_MIN:TREES_MAX, 1)
  density <- runif(1, 450, 650)                   # surviving trees per acre
  dbh     <- round(3 + rweibull(n, shape = 3.5, scale = 0.32 * age * si / 70), 1)
  ht      <- curtis_arney(dbh) * exp(rnorm(n, 0, 0.05))

  # rescale so the plot's top height (tallest 40 TPA by DBH) matches the site curve
  tpa <- rep(density / n, n)
  ord <- order(-dbh)
  w   <- pmin(tpa[ord], pmax(0, 40 - c(0, cumsum(tpa[ord])[-n])))
  top <- sum(ht[ord] * w) / sum(w)
  ht  <- 4.5 + (ht - 4.5) * (carmean(age, si) - 4.5) / (top - 4.5)
  ht  <- pmin(ht, 0.95 * CM[["c1"]] * si^CM[["c2"]])

  data.frame(PlotID = id, TreeID = seq_len(n), Dbh = dbh, Ht = round(ht, 1),
             CR = round(runif(n, 35, 60)), TPA = round(tpa, 4))
}

trees <- do.call(rbind, Map(make_trees, plots$PlotID, plots$Age, plots$SI))
rownames(trees) <- NULL

dir.create("data", showWarnings = FALSE)
write.csv(plots, "data/plots.csv", row.names = FALSE)
write.csv(trees, "data/trees.csv", row.names = FALSE)
cat(sprintf("Wrote %d plots, %d trees (%.0f per plot)\n",
            nrow(plots), nrow(trees), nrow(trees) / nrow(plots)))
