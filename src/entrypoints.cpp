#include "window_metrics.h"

using namespace Rcpp;

static SSPars make_pars(List p) {
  SSPars s;
  s.nbins   = as<int>(p["nbins"]);
  s.gmin    = as<double>(p["gmin"]);
  s.gmax    = as<double>(p["gmax"]);
  s.nlevels = as<int>(p["nlevels"]);
  s.lacbox  = as<int>(p["lacbox"]);
  s.resx    = as<double>(p["resx"]);
  s.resy    = as<double>(p["resy"]);
  return s;
}

/* gather a complete window; returns false if it touches the edge or an NA */
static inline bool gather(const NumericVector& z, int nrow, int ncol,
                          int r, int c, int hr, int hc,
                          std::vector<double>& buf) {
  if (r - hr < 0 || r + hr >= nrow || c - hc < 0 || c + hc >= ncol) return false;
  int i = 0;
  for (int rr = r - hr; rr <= r + hr; ++rr) {
    const R_xlen_t base = (R_xlen_t) rr * ncol;
    for (int cc = c - hc; cc <= c + hc; ++cc) {
      const double v = z[base + cc];
      if (ISNAN(v)) return false;
      buf[i++] = v;
    }
  }
  return true;
}

// [[Rcpp::export]]
NumericMatrix ss_focal_metrics_cpp(NumericVector z, int nrow, int ncol,
                                   int wr, int wc, IntegerVector mids,
                                   List pars) {
  const SSPars p = make_pars(pars);
  const std::vector<int> ids(mids.begin(), mids.end());
  const R_xlen_t n = (R_xlen_t) nrow * ncol;
  const int nm = ids.size();
  NumericMatrix out(n, nm);
  std::fill(out.begin(), out.end(), NA_REAL);

  const int hr = wr / 2, hc = wc / 2;
  std::vector<double> buf((size_t) wr * wc);
  std::vector<double> res(nm);

  for (int r = 0; r < nrow; ++r) {
    Rcpp::checkUserInterrupt();
    for (int c = 0; c < ncol; ++c) {
      if (!gather(z, nrow, ncol, r, c, hr, hc, buf)) continue;
      ss_window_metrics(buf.data(), wr, wc, ids, p, res.data());
      const R_xlen_t idx = (R_xlen_t) r * ncol + c;
      for (int m = 0; m < nm; ++m) out(idx, m) = res[m];
    }
  }
  return out;
}

/*
 * Within-window decomposition of change.
 *
 * For each cell the two dated windows W1 and W2 are gathered and two
 * counterfactual windows are built by rank transfer:
 *   A = arrangement of W2 with the value multiset of W1
 *   B = arrangement of W1 with the value multiset of W2
 * The Shapley-symmetrised components are then
 *   arrangement = 0.5 * [ (M(A) - M(W1)) + (M(W2) - M(B)) ]
 *   marginal    = 0.5 * [ (M(W2) - M(A)) + (M(B) - M(W1)) ]
 * and arrangement + marginal = M(W2) - M(W1) exactly.
 *
 * Columns per metric: total, marginal, arrangement.
 */
// [[Rcpp::export]]
NumericMatrix ss_decompose_cpp(NumericVector z1, NumericVector z2,
                               int nrow, int ncol, int wr, int wc,
                               IntegerVector mids, List pars) {
  const SSPars p = make_pars(pars);
  const std::vector<int> ids(mids.begin(), mids.end());
  const R_xlen_t n = (R_xlen_t) nrow * ncol;
  const int nm = ids.size();
  NumericMatrix out(n, 3 * nm);
  std::fill(out.begin(), out.end(), NA_REAL);

  const int hr = wr / 2, hc = wc / 2, k = wr * wc;
  std::vector<double> b1(k), b2(k), ba(k), bb(k), srt(k);
  std::vector<int> ord(k);
  std::vector<double> M1(nm), M2(nm), MA(nm), MB(nm);

  for (int r = 0; r < nrow; ++r) {
    Rcpp::checkUserInterrupt();
    for (int c = 0; c < ncol; ++c) {
      if (!gather(z1, nrow, ncol, r, c, hr, hc, b1)) continue;
      if (!gather(z2, nrow, ncol, r, c, hr, hc, b2)) continue;

      ss_rank_transfer(b2.data(), b1.data(), k, srt, ord, ba.data());
      ss_rank_transfer(b1.data(), b2.data(), k, srt, ord, bb.data());

      ss_window_metrics(b1.data(), wr, wc, ids, p, M1.data());
      ss_window_metrics(b2.data(), wr, wc, ids, p, M2.data());
      ss_window_metrics(ba.data(), wr, wc, ids, p, MA.data());
      ss_window_metrics(bb.data(), wr, wc, ids, p, MB.data());

      const R_xlen_t idx = (R_xlen_t) r * ncol + c;
      for (int m = 0; m < nm; ++m) {
        const double tot = M2[m] - M1[m];
        const double arr = 0.5 * ((MA[m] - M1[m]) + (M2[m] - MB[m]));
        const double mar = 0.5 * ((M2[m] - MA[m]) + (MB[m] - M1[m]));
        out(idx, m)            = tot;
        out(idx, nm + m)       = mar;
        out(idx, 2 * nm + m)   = arr;
      }
    }
  }
  return out;
}

/*
 * Spatial permutation test at sampled windows.
 *
 * For each sampled window the metrics are evaluated on the observed window
 * and on `nperm` random permutations of its values. The result has
 * (1 + nperm) * ns rows: the first ns rows hold the observed values and each
 * following block of ns rows holds one permutation. Windows that touch the
 * raster edge or contain a missing value return NA.
 */
// [[Rcpp::export]]
NumericMatrix ss_permute_cpp(NumericVector z, int nrow, int ncol,
                             int wr, int wc, IntegerVector mids,
                             IntegerVector cells, int nperm, List pars) {
  const SSPars p = make_pars(pars);
  const std::vector<int> ids(mids.begin(), mids.end());
  const int nm = ids.size(), k = wr * wc;
  const int hr = wr / 2, hc = wc / 2;
  const int ns = cells.size();
  NumericMatrix out((R_xlen_t) (1 + nperm) * ns, nm);
  std::fill(out.begin(), out.end(), NA_REAL);

  std::vector<double> buf(k), perm(k), res(nm);

  for (int s = 0; s < ns; ++s) {
    Rcpp::checkUserInterrupt();
    const R_xlen_t idx = cells[s];
    const int r = (int) (idx / ncol), c = (int) (idx % ncol);
    if (!gather(z, nrow, ncol, r, c, hr, hc, buf)) continue;

    ss_window_metrics(buf.data(), wr, wc, ids, p, res.data());
    for (int m = 0; m < nm; ++m) out(s, m) = res[m];

    for (int b = 1; b <= nperm; ++b) {
      perm = buf;
      for (int i = k - 1; i > 0; --i) {
        const int j = (int) (unif_rand() * (i + 1));
        std::swap(perm[i], perm[j]);
      }
      ss_window_metrics(perm.data(), wr, wc, ids, p, res.data());
      const R_xlen_t row = (R_xlen_t) b * ns + s;
      for (int m = 0; m < nm; ++m) out(row, m) = res[m];
    }
  }
  return out;
}

/*
 * Metrics computed over the whole raster treated as a single window, on the
 * observed arrangement (first row) and on `nperm` random permutations of all
 * cells (following rows). The raster must not contain missing values.
 */
// [[Rcpp::export]]
NumericMatrix ss_whole_cpp(NumericVector z, int nrow, int ncol,
                           IntegerVector mids, int nperm, List pars) {
  const SSPars p = make_pars(pars);
  const std::vector<int> ids(mids.begin(), mids.end());
  const int nm = ids.size();
  const R_xlen_t k = (R_xlen_t) nrow * ncol;
  NumericMatrix out(1 + nperm, nm);
  std::fill(out.begin(), out.end(), NA_REAL);

  std::vector<double> buf(z.begin(), z.end()), res(nm);
  for (R_xlen_t i = 0; i < k; ++i) if (ISNAN(buf[i])) return out;

  ss_window_metrics(buf.data(), nrow, ncol, ids, p, res.data());
  for (int m = 0; m < nm; ++m) out(0, m) = res[m];

  std::vector<double> perm(k);
  for (int b = 1; b <= nperm; ++b) {
    Rcpp::checkUserInterrupt();
    perm = buf;
    for (R_xlen_t i = k - 1; i > 0; --i) {
      const R_xlen_t j = (R_xlen_t) (unif_rand() * (i + 1));
      std::swap(perm[i], perm[j]);
    }
    ss_window_metrics(perm.data(), nrow, ncol, ids, p, res.data());
    for (int m = 0; m < nm; ++m) out(b, m) = res[m];
  }
  return out;
}
