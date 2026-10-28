# Project every plot forward and build one summary table.
#
# Both functions give the same table. `grow` is the C++ growth function to
# use for each plot -- grow_plot (src/grow.cpp) or grow_plot_fast
# (src/grow_fast.cpp) -- so R-side and C++-side fixes can be timed separately.

# The version the team wrote first.
project_stands_slow <- function(trees, plots, years = 30, grow = grow_plot) {
  out <- NULL
  for (k in seq_len(nrow(plots))) {
    id <- plots$PlotID[k]
    pt <- trees[trees$PlotID == id, ]
    res <- grow(pt$Dbh, pt$Ht, pt$CR, pt$TPA, plots$Age[k], plots$SI[k], years)
    out <- rbind(out, data.frame(PlotID = id, res))
  }
  out
}

# Same result: split the tree table once, collect pieces in a list, bind once.
project_stands_fast <- function(trees, plots, years = 30, grow = grow_plot) {
  rows <- split(seq_len(nrow(trees)), trees$PlotID)
  out <- vector("list", nrow(plots))
  for (k in seq_len(nrow(plots))) {
    id <- plots$PlotID[k]
    i <- rows[[id]]
    res <- grow(trees$Dbh[i], trees$Ht[i], trees$CR[i], trees$TPA[i],
                plots$Age[k], plots$SI[k], years)
    out[[k]] <- data.frame(PlotID = id, res)
  }
  do.call(rbind, out)
}
