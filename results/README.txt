Pre-generated callgrind output, produced by the "Callgrind profile" workflow
on ubuntu-latest (valgrind 3.22, R 4.6.1), 200 plots grown 30 years.

  callgrind.out.naive / .fast            raw callgrind data (open in QCacheGrind/KCachegrind)
  annotate_inclusive_naive.txt / _fast   per-function report (--inclusive=yes)
  annotate_naive.txt / _fast             per-line report (--auto=yes)

Use these if your own run failed or is queued:
  git fetch && git checkout origin/results -- results/
