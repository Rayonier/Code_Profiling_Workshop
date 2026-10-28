// Same projection as grow.cpp -- same published FVS Southern variant
// equations (fvs_sn_lp.h), same answers -- with the costly habits removed:
//   1. PBAL from one sort + running sum, instead of scanning every tree
//      for every tree (FVS itself does it this way, via its BA percentiles).
//   2. Terms that don't change from tree to tree (or, for crown ratio, from
//      cycle to cycle) are computed once.
//   3. No per-cycle copies of the tree vectors: growth is computed first,
//      then applied in place.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include "fvs_sn_lp.h"

#ifdef USE_CALLGRIND
#include <valgrind/callgrind.h>
#else
#define CALLGRIND_START_INSTRUMENTATION
#define CALLGRIND_STOP_INSTRUMENTATION
#endif

using namespace Rcpp;

//' Faster grow_plot(): identical inputs and output.
// [[Rcpp::export]]
DataFrame grow_plot_fast(NumericVector dbh_in, NumericVector ht_in, NumericVector cr,
                         NumericVector tpa_in, double age, double si, int years,
                         double slope = 0.05, double aspect = 0.0) {
  const int n = dbh_in.size();
  const int ncycles = years / static_cast<int>(lp::CYCLE);
  std::vector<double> dbh(dbh_in.begin(), dbh_in.end());
  std::vector<double> ht(ht_in.begin(), ht_in.end());
  std::vector<double> tpa(tpa_in.begin(), tpa_in.end());

  // --- terms that are the same for every tree in the plot -------------------
  // {4.7.1.1}: everything except the tree- and density-level terms
  const double ln_dds_const = lp::B1 + lp::B6 * si + lp::B9 * slope +
                              lp::B10 * std::cos(aspect) * slope +
                              lp::B11 * std::sin(aspect) * slope +
                              lp::FORTYPE + lp::ECOUNIT + lp::PLANT;
  // {4.6.1.1}: site-index parts of the Carmean curve
  const double htmax = lp::C1 * std::pow(si, lp::C2);
  const double shape = lp::C4 * std::pow(si, lp::C5);
  // {4.7.2.3}-{4.7.2.6}: shade-tolerance constants
  const double fctrkx = std::pow(lp::RHK / lp::RHYXS, lp::RHM - 1.0) - 1.0;
  const double fctrrb = -lp::RHR / (1.0 - lp::RHB);
  const double rhxs_term = std::pow(lp::RHXS, 1.0 - lp::RHB);
  const double fctrm = 1.0 / (1.0 - lp::RHM);

  // crown ratio never changes, so its terms are computed once per tree
  std::vector<double> ln_cr(n), hgmdcr(n);
  for (int i = 0; i < n; ++i) {
    const double crp = cr[i] / 100.0;
    ln_cr[i] = std::log(cr[i]);
    hgmdcr[i] = std::min(100.0 * crp * crp * crp * std::exp(-5.0 * crp), 1.0);
  }

  std::vector<int> ord(n);
  std::vector<double> pbal(n), dg(n), htg(n), surv(n);
  const int nrow = ncycles + 1;
  NumericVector out_year(nrow), out_age(nrow), out_tpa(nrow), out_ba(nrow),
                out_qmd(nrow), out_hd(nrow);

  CALLGRIND_START_INSTRUMENTATION;
  for (int cycle = 0; cycle <= ncycles; ++cycle) {
    // one sort by DBH (largest first) gives top height and PBAL
    std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(),
                     [&](int a, int b) { return dbh[a] > dbh[b]; });

    double sum_tpa = 0.0, sum_ht = 0.0;
    for (int i : ord) {
      double p = tpa[i];
      if (sum_tpa + p > lp::TOP_TPA) p = lp::TOP_TPA - sum_tpa;
      sum_tpa += p;
      sum_ht += ht[i] * p;
      if (sum_tpa >= lp::TOP_TPA) break;
    }
    const double avh = sum_tpa > 0.0 ? sum_ht / sum_tpa : 0.0;

    // running BA total down the sorted list; trees tied on DBH share the
    // total of the strictly larger trees before them
    double ba = 0.0, trees = 0.0;
    for (int k = 0; k < n; ) {
      int k_end = k;
      while (k_end < n && dbh[ord[k_end]] == dbh[ord[k]]) ++k_end;
      double group_ba = 0.0;
      for (int m = k; m < k_end; ++m) {
        int i = ord[m];
        pbal[i] = ba;
        group_ba += lp::BA_FACTOR * dbh[i] * dbh[i] * tpa[i];
        trees += tpa[i];
      }
      ba += group_ba;
      k = k_end;
    }

    out_year[cycle] = cycle * lp::CYCLE;
    out_age[cycle] = age + cycle * lp::CYCLE;
    out_tpa[cycle] = trees;
    out_ba[cycle] = ba;
    out_qmd[cycle] = std::sqrt(ba / (trees * lp::BA_FACTOR));
    out_hd[cycle] = avh;
    if (cycle == ncycles) break;

    const double ln_dds_stand = ln_dds_const + lp::B7 * ba;
    for (int i = 0; i < n; ++i) {
      const double d = dbh[i], h = ht[i];
      const double relht = avh > 0.0 ? std::min(h / avh, lp::RELHT_MAX) : 0.0;

      // {4.7.1.1} diameter growth
      double ln_dds = ln_dds_stand + lp::B2 * std::log(d) + lp::B3 * d * d +
                      lp::B4 * ln_cr[i] + lp::B5 * relht + lp::B8 * pbal[i];
      ln_dds = std::max(ln_dds, lp::LN_DDS_MIN);
      const double bratio = lp::bark_ratio(d);
      const double dib = d * bratio;
      dg[i] = (std::sqrt(dib * dib + std::exp(ln_dds)) - dib) / bratio;

      // {4.6.1.1}, {4.7.2.1}-{4.7.2.7} height growth
      double potential = 0.0;
      if (htmax - h > 1.0) {
        const double aget = std::log(1.0 - std::pow(h / htmax, 1.0 / shape)) / lp::C3;
        potential = htmax * std::pow(1.0 - std::exp(lp::C3 * (aget + lp::CYCLE)), shape) - h;
      }
      const double fctrxb = std::pow(relht, 1.0 - lp::RHB) - rhxs_term;
      const double hgmdrh = lp::RHK * std::pow(1.0 + fctrkx * std::exp(fctrrb * fctrxb), fctrm);
      htg[i] = std::max(potential * (lp::WTCR * hgmdcr[i] + (1.0 - lp::WTCR) * hgmdrh), lp::HTG_MIN);

      // {5.0.1}-{5.0.2} background mortality; survival over 5 years is
      // (1 - RI)^5, multiplied out rather than calling pow()
      const double ri = 1.0 / (1.0 + std::exp(lp::P0 + lp::P1 * d));
      const double s = 1.0 - ri;
      surv[i] = s * s * s * s * s;
    }
    for (int i = 0; i < n; ++i) {
      dbh[i] += dg[i];
      ht[i] += htg[i];
      tpa[i] *= surv[i];
    }
  }
  CALLGRIND_STOP_INSTRUMENTATION;

  return DataFrame::create(_["Year"] = out_year, _["Age"] = out_age,
                           _["TPA"] = out_tpa, _["BA"] = out_ba,
                           _["QMD"] = out_qmd, _["HD"] = out_hd);
}
