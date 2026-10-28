// Individual-tree growth projection for a loblolly pine plot, using the
// published FVS Southern variant equations in fvs_sn_lp.h.
//
// This is the version the team wrote first. It gives correct answers --
// but it has a few habits that cost time. Section 4 of the workshop uses
// callgrind to find them.

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

// Tree order by DBH, largest first (ties keep input order).
static std::vector<int> order_by_dbh(const NumericVector& dbh) {
  std::vector<int> ord(dbh.size());
  std::iota(ord.begin(), ord.end(), 0);
  std::stable_sort(ord.begin(), ord.end(),
                   [&](int a, int b) { return dbh[a] > dbh[b]; });
  return ord;
}

// FVS top height (AVH): mean height of the largest trees making up 40 TPA.
static double top_height(const NumericVector& dbh, const NumericVector& ht,
                         const NumericVector& tpa) {
  std::vector<int> ord = order_by_dbh(dbh);
  double sum_tpa = 0.0, sum_ht = 0.0;
  for (int i : ord) {
    double p = tpa[i];
    if (sum_tpa + p > lp::TOP_TPA) p = lp::TOP_TPA - sum_tpa;
    sum_tpa += p;
    sum_ht += ht[i] * p;
    if (sum_tpa >= lp::TOP_TPA) break;
  }
  return sum_tpa > 0.0 ? sum_ht / sum_tpa : 0.0;
}

// {4.7.1.1} 5-year diameter growth (outside bark, inches).
static double diameter_growth(double dbh, double ht, double cr, double avh,
                              double si, double ba, double pbal,
                              double slope, double aspect) {
  double relht = avh > 0.0 ? std::min(ht / avh, lp::RELHT_MAX) : 0.0;
  double ln_dds = lp::B1 + lp::B2 * std::log(dbh) + lp::B3 * dbh * dbh +
                  lp::B4 * std::log(cr) + lp::B5 * relht + lp::B6 * si +
                  lp::B7 * ba + lp::B8 * pbal + lp::B9 * slope +
                  lp::B10 * std::cos(aspect) * slope +
                  lp::B11 * std::sin(aspect) * slope +
                  lp::FORTYPE + lp::ECOUNIT + lp::PLANT;
  ln_dds = std::max(ln_dds, lp::LN_DDS_MIN);
  double bratio = lp::bark_ratio(dbh);
  double dib = dbh * bratio;                       // inside-bark diameter
  return (std::sqrt(dib * dib + std::exp(ln_dds)) - dib) / bratio;
}

// {4.6.1.1} and {4.7.2.1}-{4.7.2.7} 5-year height growth (feet).
static double height_growth(double ht, double cr, double avh, double si) {
  double htmax = lp::C1 * std::pow(si, lp::C2);
  double potential = 0.0;
  if (htmax - ht > 1.0) {
    double shape = lp::C4 * std::pow(si, lp::C5);
    double aget = std::log(1.0 - std::pow(ht / htmax, 1.0 / shape)) / lp::C3;
    double ht5 = htmax * std::pow(1.0 - std::exp(lp::C3 * (aget + lp::CYCLE)), shape);
    potential = ht5 - ht;
  }
  double crp = cr / 100.0;
  double hgmdcr = std::min(100.0 * crp * crp * crp * std::exp(-5.0 * crp), 1.0);
  double relht = avh > 0.0 ? std::min(ht / avh, lp::RELHT_MAX) : 0.0;
  double fctrkx = std::pow(lp::RHK / lp::RHYXS, lp::RHM - 1.0) - 1.0;
  double fctrrb = -lp::RHR / (1.0 - lp::RHB);
  double fctrxb = std::pow(relht, 1.0 - lp::RHB) - std::pow(lp::RHXS, 1.0 - lp::RHB);
  double fctrm = 1.0 / (1.0 - lp::RHM);
  double hgmdrh = lp::RHK * std::pow(1.0 + fctrkx * std::exp(fctrrb * fctrxb), fctrm);
  double htg = potential * (lp::WTCR * hgmdcr + (1.0 - lp::WTCR) * hgmdrh);
  return std::max(htg, lp::HTG_MIN);
}

// {5.0.1}-{5.0.2} background mortality rate over one cycle.
static double mortality_rate(double dbh) {
  double ri = 1.0 / (1.0 + std::exp(lp::P0 + lp::P1 * dbh));
  return 1.0 - std::pow(1.0 - ri, lp::CYCLE);
}

//' Project one plot's tree list forward in 5-year cycles.
//' Returns one row per cycle (including the start): Year, Age, TPA, BA, QMD, HD.
// [[Rcpp::export]]
DataFrame grow_plot(NumericVector dbh, NumericVector ht, NumericVector cr,
                    NumericVector tpa, double age, double si, int years,
                    double slope = 0.05, double aspect = 0.0) {
  int n = dbh.size();
  int ncycles = years / static_cast<int>(lp::CYCLE);
  dbh = clone(dbh); ht = clone(ht); tpa = clone(tpa);    // don't modify the caller's data

  std::vector<double> out_year, out_age, out_tpa, out_ba, out_qmd, out_hd;

  CALLGRIND_START_INSTRUMENTATION;
  for (int cycle = 0; cycle <= ncycles; ++cycle) {
    // keep this cycle's starting state; every tree grows from the same snapshot
    NumericVector dbh0 = clone(dbh), ht0 = clone(ht), tpa0 = clone(tpa);

    double ba = 0.0, trees = 0.0;
    for (int i = 0; i < n; ++i) {
      ba += lp::BA_FACTOR * dbh0[i] * dbh0[i] * tpa0[i];
      trees += tpa0[i];
    }
    double avh = top_height(dbh0, ht0, tpa0);

    out_year.push_back(cycle * lp::CYCLE);
    out_age.push_back(age + cycle * lp::CYCLE);
    out_tpa.push_back(trees);
    out_ba.push_back(ba);
    out_qmd.push_back(std::sqrt(ba / (trees * lp::BA_FACTOR)));
    out_hd.push_back(avh);
    if (cycle == ncycles) break;

    for (int i = 0; i < n; ++i) {
      // plot basal area in larger trees: add up every tree bigger than this one
      double pbal = 0.0;
      for (int j = 0; j < n; ++j) {
        if (dbh0[j] > dbh0[i]) pbal += lp::BA_FACTOR * dbh0[j] * dbh0[j] * tpa0[j];
      }
      dbh[i] = dbh0[i] + diameter_growth(dbh0[i], ht0[i], cr[i], avh, si, ba, pbal,
                                         slope, aspect);
      ht[i] = ht0[i] + height_growth(ht0[i], cr[i], avh, si);
      tpa[i] = tpa0[i] * (1.0 - mortality_rate(dbh0[i]));
    }
  }
  CALLGRIND_STOP_INSTRUMENTATION;

  return DataFrame::create(_["Year"] = out_year, _["Age"] = out_age,
                           _["TPA"] = out_tpa, _["BA"] = out_ba,
                           _["QMD"] = out_qmd, _["HD"] = out_hd);
}
