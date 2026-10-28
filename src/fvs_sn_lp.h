// Published equations and coefficients for loblolly pine (FVS species code LP)
// from the Southern (SN) variant of the Forest Vegetation Simulator.
//
// Source: FVS Staff. 2008 (revised March 30, 2026). Southern (SN) Variant
//   Overview - Forest Vegetation Simulator. USDA Forest Service, Forest
//   Management Service Center, Fort Collins, CO.
//   https://www.fs.usda.gov/fmsc/ftp/fvs/docs/overviews/FVSsn_Overview.pdf
// Coefficients cross-checked against the open-source FVS code
// (github.com/USDAForestService/ForestVegetationSimulator, sn/*.f).
// {x.y.z} numbers below refer to equations in the overview.
//
// This is a teaching re-implementation of selected equations, NOT FVS.
// Simplifications: crown ratio held constant; large-tree models only
// (DBH >= 3 in); background mortality only (no SDI-based density mortality
// or mortality dispersal); no random error, calibration or regeneration.

#ifndef FVS_SN_LP_H
#define FVS_SN_LP_H

namespace lp {

const double CYCLE = 5.0;               // FVS projection cycle length (years)
const double BA_FACTOR = 0.005454154;   // ft^2 of basal area per in^2 of DBH^2

// {4.7.1.1} large-tree diameter growth, Tables 4.7.1.1 - 4.7.1.7
const double B1 = 0.222214, B2 = 1.163040, B3 = -0.000863, B4 = 0.028483;
const double B5 = 0.006935, B6 = 0.005018, B7 = -0.003408, B8 = -0.004184;
const double B9 = -0.759347, B10 = 0.185360, B11 = -0.072842;
const double FORTYPE = 0.0;             // base forest type group FTYLPN
const double ECOUNIT = 0.0;             // base ecological unit group 232
const double PLANT = 0.245669;          // managed pine plantation
const double LN_DDS_MIN = -9.21;        // floor applied in FVS dgf.f
// {4.7.1.3} bounding function: LP limits are 998/999 in, so it never applies.

// {4.2.1} bark ratio, Table 4.2.1
const double BARK1 = -0.48140, BARK2 = 0.91413;

// {4.6.1.1} potential height (Carmean et al. 1989), Table 4.6.1.1
const double C1 = 1.1421, C2 = 1.0042, C3 = -0.0374, C4 = 0.7632, C5 = 0.0358;

// {4.7.2.1}-{4.7.2.7} large-tree height growth modifiers, Tables 4.7.2.1-2
// (LP is "Intolerant")
const double WTCR = 0.25;               // weight on crown ratio modifier
const double RHR = 13.0, RHYXS = 0.05, RHM = 1.1, RHB = -1.60, RHXS = 0.0, RHK = 1.0;
const double RELHT_MAX = 1.5;           // cap applied in FVS dgf.f / htgf.f
const double HTG_MIN = 0.1;             // floor applied in FVS htgf.f

// {5.0.1}-{5.0.2} background mortality, Table 5.0.1
const double P0 = 5.5876999, P1 = -0.0053480;

// FVS top height (AVH): TPA-weighted mean height of the largest-diameter
// trees making up this many trees per acre (FVS avht40.f).
const double TOP_TPA = 40.0;

inline double bark_ratio(double dob) {
  double r = (BARK1 + BARK2 * dob) / dob;
  if (r > 0.99) r = 0.99;
  if (r < 0.80) r = 0.80;
  return r;
}

}  // namespace lp

#endif
