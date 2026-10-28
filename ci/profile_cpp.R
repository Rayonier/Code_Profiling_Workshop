# Driver that the GitHub workflow runs under valgrind/callgrind.
# Callgrind only records between the CALLGRIND_START/STOP markers inside the
# C++, so R's own startup, the compile, and reading the CSV aren't counted.
#
# Settings (environment variables):
#   CPP_VERSION  naive (src/grow.cpp, default) or fast (src/grow_fast.cpp)
#   N_PLOTS      plots to grow (default 200; callgrind runs ~50x slower than normal)
#
# Runs fine without valgrind too:  Rscript ci/profile_cpp.R

version <- Sys.getenv("CPP_VERSION", "naive")
n_plots <- as.integer(Sys.getenv("N_PLOTS", "200"))
src <- switch(version, naive = "src/grow.cpp", fast = "src/grow_fast.cpp",
              stop("CPP_VERSION must be 'naive' or 'fast', not '", version, "'"))

Rcpp::sourceCpp(src)
grow <- if (version == "fast") grow_plot_fast else grow_plot

plots <- read.csv("data/plots.csv")[seq_len(n_plots), ]
# trees.csv is sorted by plot, so the first plots' trees are at the top of the file
trees <- read.csv("data/trees.csv", nrows = n_plots * 260)
trees <- trees[trees$PlotID %in% plots$PlotID, ]
stopifnot(setequal(unique(trees$PlotID), plots$PlotID))

rows <- split(seq_len(nrow(trees)), trees$PlotID)
t0 <- proc.time()[["elapsed"]]
for (k in seq_len(nrow(plots))) {
  i <- rows[[plots$PlotID[k]]]
  grow(trees$Dbh[i], trees$Ht[i], trees$CR[i], trees$TPA[i],
       plots$Age[k], plots$SI[k], 30)
}
cat(sprintf("Grew %d plots (%d trees) with the %s C++ in %.1f s\n",
            nrow(plots), nrow(trees), version, proc.time()[["elapsed"]] - t0))
