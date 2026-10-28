# Maintainer checks: model correctness, naive vs fast equality, and the timing
# story the workshop depends on. Run from the workshop folder:
#   Rscript dev/check.R          (all checks, ~1 minute)
#   Rscript dev/check.R quick    (skip the full-size timings)

args <- commandArgs(trailingOnly = TRUE)
quick <- "quick" %in% args

Rcpp::sourceCpp("src/grow.cpp")
Rcpp::sourceCpp("src/grow_fast.cpp")
source("R/project_stands.R")
trees <- read.csv("data/trees.csv")
plots <- read.csv("data/plots.csv")

ok <- function(cond, msg) {
  cat(if (isTRUE(cond)) "PASS " else "FAIL ", msg, "\n", sep = "")
  if (!isTRUE(cond)) quit(status = 1)
}

# --- 1. Independent R implementation of the same published equations --------
# Written separately from the C++ (vectorized, PBAL via outer()) so a mistake
# in one is unlikely to be repeated in the other. Coefficients: FVS-SN LP.
reference_plot <- function(dbh, ht, cr, tpa, age, si, years, slope = 0.05, aspect = 0) {
  b <- c(0.222214, 1.163040, -0.000863, 0.028483, 0.006935, 0.005018,
         -0.003408, -0.004184, -0.759347, 0.185360, -0.072842)
  plant <- 0.245669
  cc <- c(1.1421, 1.0042, -0.0374, 0.7632, 0.0358)
  baf <- 0.005454154
  top_height <- function(dbh, ht, tpa) {
    o <- order(-dbh)                          # order() is stable, like std::stable_sort
    w <- pmin(tpa[o], pmax(0, 40 - c(0, cumsum(tpa[o])[-length(o)])))
    sum(ht[o] * w) / sum(w)
  }
  rows <- list()
  for (cycle in 0:(years / 5)) {
    ba <- sum(baf * dbh^2 * tpa)
    avh <- top_height(dbh, ht, tpa)
    rows[[cycle + 1]] <- c(Year = cycle * 5, Age = age + cycle * 5, TPA = sum(tpa),
                           BA = ba, QMD = sqrt(ba / (sum(tpa) * baf)), HD = avh)
    if (cycle == years / 5) break
    pbal <- colSums(outer(dbh, dbh, ">") * (baf * dbh^2 * tpa))  # [j, i]: tree j larger than i
    relht <- pmin(ht / avh, 1.5)
    ln_dds <- b[1] + b[2] * log(dbh) + b[3] * dbh^2 + b[4] * log(cr) + b[5] * relht +
      b[6] * si + b[7] * ba + b[8] * pbal + b[9] * slope + b[10] * cos(aspect) * slope +
      b[11] * sin(aspect) * slope + plant
    ln_dds <- pmax(ln_dds, -9.21)
    br <- pmin(pmax((-0.48140 + 0.91413 * dbh) / dbh, 0.80), 0.99)
    dib <- dbh * br
    dg <- (sqrt(dib^2 + exp(ln_dds)) - dib) / br
    htmax <- cc[1] * si^cc[2]; shape <- cc[4] * si^cc[5]
    pot <- ifelse(htmax - ht > 1, {
      aget <- log(1 - (pmin(ht, htmax) / htmax)^(1 / shape)) / cc[3]
      htmax * (1 - exp(cc[3] * (aget + 5)))^shape - ht
    }, 0)
    crp <- cr / 100
    hgmdcr <- pmin(100 * crp^3 * exp(-5 * crp), 1)
    hgmdrh <- (1 + ((1 / 0.05)^0.1 - 1) * exp((-13 / 2.6) * relht^2.6))^(1 / (1 - 1.1))
    htg <- pmax(pot * (0.25 * hgmdcr + 0.75 * hgmdrh), 0.1)
    ri <- 1 / (1 + exp(5.5876999 - 0.0053480 * dbh))
    tpa <- tpa * (1 - ri)^5
    dbh <- dbh + dg
    ht <- ht + htg
  }
  as.data.frame(do.call(rbind, rows))
}

set.seed(1)
for (k in sample(nrow(plots), 5)) {
  pt <- trees[trees$PlotID == plots$PlotID[k], ]
  args_k <- list(pt$Dbh, pt$Ht, pt$CR, pt$TPA, plots$Age[k], plots$SI[k], 30)
  ref <- do.call(reference_plot, args_k)
  ok(isTRUE(all.equal(do.call(grow_plot, args_k), ref, tolerance = 1e-8, check.attributes = FALSE)),
     sprintf("grow_plot matches the R reference implementation (plot %s)", plots$PlotID[k]))
}

# --- 2. Model sanity --------------------------------------------------------
carmean <- function(age, si) 1.1421 * si^1.0042 * (1 - exp(-0.0374 * age))^(0.7632 * si^0.0358)
res <- project_stands_fast(trees, plots, years = 30, grow = grow_plot_fast)
res <- merge(res, plots[, c("PlotID", "SI")], by = "PlotID")
dev <- res$HD - carmean(res$Age, res$SI)
ok(all(is.finite(as.matrix(res[, c("TPA", "BA", "QMD", "HD")]))), "no NaN/Inf in any projection")
ok(max(abs(dev[res$Year == 0])) < 0.5, "starting top height matches the Carmean site curve")
cat(sprintf("     top height vs site curve after 30 yr: mean %+.1f ft, range %+.1f to %+.1f ft\n",
            mean(dev[res$Year == 30]), min(dev[res$Year == 30]), max(dev[res$Year == 30])))
ok(abs(mean(dev[res$Year == 30])) < 5, "projected top height stays within 5 ft of the site curve on average")
end <- res[res$Year == 30, ]
cat(sprintf("     at +30 yr: QMD %.1f-%.1f in, BA %.0f-%.0f ft2/ac, TPA %.0f-%.0f\n",
            min(end$QMD), max(end$QMD), min(end$BA), max(end$BA), min(end$TPA), max(end$TPA)))
ok(all(diff(res$QMD[order(res$PlotID, res$Year)])[res$Year[order(res$PlotID, res$Year)][-1] > 0] > 0),
   "QMD increases every cycle")

# near the height asymptote: no NaN, minimum growth applies
htmax <- 1.1421 * 60^1.0042
edge <- grow_plot(c(8, 9), c(htmax - 0.5, htmax - 0.1), c(40, 40), c(300, 300), 60, 60, 10)
ok(all(is.finite(as.matrix(edge))), "trees at the height asymptote don't produce NaN")

# --- 3. Naive vs fast: same answers ------------------------------------------
sub <- plots[seq_len(if (quick) 100 else nrow(plots)), ]
a <- project_stands_slow(trees, sub)
b <- project_stands_fast(trees, sub)
f <- project_stands_fast(trees, sub, grow = grow_plot_fast)
ok(isTRUE(all.equal(a, b)), "project_stands_slow == project_stands_fast")
ok(isTRUE(all.equal(b, f, tolerance = 1e-10)), "grow_plot == grow_plot_fast (1e-10)")

# --- 4. The timing story -----------------------------------------------------
if (!quick) {
  rows <- split(seq_len(nrow(trees)), trees$PlotID)
  cpp_only <- function(g) for (k in seq_len(nrow(plots))) {
    i <- rows[[plots$PlotID[k]]]
    g(trees$Dbh[i], trees$Ht[i], trees$CR[i], trees$TPA[i], plots$Age[k], plots$SI[k], 30)
  }
  el <- function(expr) system.time(expr)[["elapsed"]]
  t <- c(slow_R_naive_cpp = el(project_stands_slow(trees, plots)),
         fast_R_naive_cpp = el(project_stands_fast(trees, plots)),
         fast_R_fast_cpp  = el(project_stands_fast(trees, plots, grow = grow_plot_fast)),
         naive_cpp_only   = el(cpp_only(grow_plot)),
         fast_cpp_only    = el(cpp_only(grow_plot_fast)))
  print(round(t, 2))
  share_slow  <- t[["naive_cpp_only"]] / t[["slow_R_naive_cpp"]]
  share_fixed <- t[["naive_cpp_only"]] / t[["fast_R_naive_cpp"]]
  cat(sprintf("     C++ share: %.0f%% of the slow pipeline, %.0f%% after the R fix\n",
              100 * share_slow, 100 * share_fixed))
  ok(t[["slow_R_naive_cpp"]] < 60, "slow pipeline finishes in under a minute")
  ok(share_slow < 0.35, "in the slow pipeline, R (not C++) is most of the time")
  ok(share_fixed > 0.5, "after the R fix, C++ is most of the time")
  ok(t[["fast_cpp_only"]] < t[["naive_cpp_only"]] * 0.75, "the C++ fix is a real speedup")
}
cat("All checks passed.\n")
