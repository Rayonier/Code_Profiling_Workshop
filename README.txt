Measure Before You Touch It: Profiling R and Rcpp Code
======================================================

Training repository for a hands-on workshop on profiling R and Rcpp code.

Workshop materials (provided by the instructor, not stored in this repo):
1. intro_deck.pptx   - Part 1 lecture slides
2. workshop.html     - Part 2 hands-on workshop (open in a browser)
3. takeaways.html    - one-page summary

Before the session: make sure you can open Rayonier/profiling-workshop,
create a branch there, and see its Actions tab. In section 1 you'll create
your own branch off main and run the "Callgrind profile" workflow on it.
Please don't push to main or results.

To run the workshop code yourself, open workshop.Rmd in RStudio with this
folder as the working directory and click Knit (or run chunks one by one).
R packages needed: Rcpp, profvis, bench, <dplyr|data.table>, rmarkdown.
Compiling src/*.cpp needs Rtools on Windows.

The callgrind step runs on GitHub's Linux runners (valgrind does not run on
Windows). If your run fails or queues, the repo's `results` branch has
pre-generated output in results/.

Plan B -- GitHub not reachable at all: ask the instructor for the offline
bundle. It contains the full workshop plus results/, and every section can
be completed from it without GitHub.

Data (simulated -- regenerate with data/make_data.R, fixed seed):
  data/plots.csv   <n> plots: PlotID, Age, HD (ft), BA (sq ft/ac), TPA, ...
