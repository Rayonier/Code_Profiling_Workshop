Measure Before You Touch It: Profiling R and Rcpp Code
======================================================

Training repository for a hands-on workshop on profiling R and Rcpp code.

Workshop materials (provided by the instructor, not stored in this repo):
1. intro_deck.pptx   - Part 1 lecture slides
2. workshop.html     - Part 2 hands-on workshop (open in a browser)
3. takeaways.html    - one-page summary
workshop.Rmd and takeaways.Rmd in this repo are their sources.

Setup (once): open R in this folder and run
    source("bootstrap.R")
It installs any missing packages (Rcpp, profvis, bench, rmarkdown, knitr,
ggplot2, pkgbuild) and checks that C++ compiles. On Windows that needs
Rtools: https://cran.r-project.org/bin/windows/Rtools/

Before the session: make sure you can open Rayonier/profiling-workshop,
create a branch there, and see its Actions tab. In section 1 you'll create
your own branch off main and run the "Callgrind profile" workflow on it.
Please don't push to main or results.

To run the workshop code yourself, open workshop.Rmd in RStudio with this
folder as the working directory and click Knit (or run chunks one by one).
Knitting takes a couple of minutes: it really runs the slow version.

The callgrind step runs on GitHub's Linux runners (valgrind does not run on
Windows). If your run fails or queues, the repo's `results` branch has
pre-generated output in results/.

Plan B -- GitHub not reachable at all: ask the instructor for the offline
bundle. It contains the full workshop plus results/, and every section can
be completed from it without GitHub.

Files:
  R/project_stands.R   the pipeline: project_stands_slow() and _fast()
  src/grow.cpp         C++ growth model as first written (grow_plot)
  src/grow_fast.cpp    same model, optimized (grow_plot_fast)
  src/fvs_sn_lp.h      published equations' coefficients, with sources
  ci/                  what the GitHub workflow runs under callgrind
  dev/check.R          maintainer checks (model, naive == fast, timings)

Data (simulated -- regenerate with `Rscript data/make_data.R`, fixed seed):
  data/plots.csv   PlotID, Age (years), SI (site index, ft at base age 50)
  data/trees.csv   PlotID, TreeID, Dbh (in), Ht (ft), CR (crown ratio, %),
                   TPA (trees per acre each record represents)

Model source:
  The growth model uses published equations and coefficients for loblolly
  pine from FVS Staff (2008, revised March 30, 2026), Southern (SN) Variant
  Overview - Forest Vegetation Simulator, USDA Forest Service:
  https://www.fs.usda.gov/fmsc/ftp/fvs/docs/overviews/FVSsn_Overview.pdf
  It is a simplified teaching re-implementation, not FVS: crown ratio is
  held constant, only the large-tree models are used, mortality is
  background mortality only (so basal area runs high late in a projection),
  and there is no random error, calibration or regeneration. The plots are
  simulated; only their tree heights come from published equations.
