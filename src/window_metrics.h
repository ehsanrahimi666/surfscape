#ifndef SURFSCAPE_WINDOW_METRICS_H
#define SURFSCAPE_WINDOW_METRICS_H

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

/*
 * Shared window-metric kernel.
 *
 * Every metric in surfscape is computed from a single w x w window held as a
 * dense row-major buffer. The same kernel is reused by three entry points:
 *   - ss_focal_metrics_cpp   : metric maps over a raster
 *   - ss_decompose_cpp       : within-window rank transfer and decomposition
 *   - ss_sensitivity_cpp     : permutation audit of configurational sensitivity
 *
 * That reuse is not an optimisation, it is a requirement of the method: the
 * decomposition evaluates metrics on counterfactual windows that never exist
 * as a raster, so they cannot be passed to an external metric package.
 *
 * Metric ids
 *   composition family (functions of the window value multiset)
 *     1 smean  2 sa    3 sq    4 ssk   5 sku   6 sp    7 sv
 *     8 sz     9 s10z 10 scv  11 sent 12 srao
 *   configuration family (depend on the arrangement of values in the window)
 *    20 sr1   21 sdq  22 sdr  23 slac
 *    30 gcon  31 gdis 32 ghom 33 gasm 34 gent  35 gcor
 */

struct SSPars {
  int    nbins   = 32;    // bins for sent / srao
  double gmin    = 0.0;   // global value range, keeps windows comparable
  double gmax    = 1.0;
  int    nlevels = 16;    // grey levels for the GLCM
  int    lacbox  = 3;     // gliding box side for slac
  double resx    = 1.0;
  double resy    = 1.0;
};

inline bool ss_needs(const std::vector<int>& ids, int lo, int hi) {
  for (size_t i = 0; i < ids.size(); ++i)
    if (ids[i] >= lo && ids[i] <= hi) return true;
  return false;
}

/*
 * buf : row-major w x w window, must be complete (no NA)
 * out : one slot per requested id
 */
inline void ss_window_metrics(const double* buf, int wr, int wc,
                              const std::vector<int>& ids,
                              const SSPars& p,
                              double* out) {

  const int k = wr * wc;
  const bool need_glcm = ss_needs(ids, 30, 35);
  const bool need_bin  = ss_needs(ids, 11, 12);

  double s1 = 0.0, mn = buf[0], mx = buf[0];
  for (int i = 0; i < k; ++i) {
    s1 += buf[i];
    if (buf[i] < mn) mn = buf[i];
    if (buf[i] > mx) mx = buf[i];
  }
  const double mean = s1 / k;

  double m2 = 0.0, m3 = 0.0, m4 = 0.0, sad = 0.0;
  for (int i = 0; i < k; ++i) {
    const double d = buf[i] - mean, d2 = d * d;
    m2 += d2; m3 += d2 * d; m4 += d2 * d2; sad += std::fabs(d);
  }
  m2 /= k; m3 /= k; m4 /= k; sad /= k;
  const double sq = std::sqrt(m2);

  /* ---- order statistics -------------------------------------------- */
  double s10z = NA_REAL;
  if (ss_needs(ids, 9, 9) && k >= 10) {
    std::vector<double> t(buf, buf + k);
    std::nth_element(t.begin(), t.begin() + 5, t.end());
    double lo = 0.0; for (int i = 0; i < 5; ++i) lo += t[i];
    std::nth_element(t.begin(), t.begin() + (k - 5), t.end());
    double hi = 0.0; for (int i = k - 5; i < k; ++i) hi += t[i];
    s10z = (hi - lo) / 5.0;
  }

  /* ---- binned distribution: Shannon entropy and Rao's Q ------------- */
  double ent = NA_REAL, rao = NA_REAL;
  if (need_bin) {
    const int nb = p.nbins;
    const double bw = (p.gmax > p.gmin) ? (p.gmax - p.gmin) / nb : 1.0;
    std::vector<double> h(nb, 0.0);
    for (int i = 0; i < k; ++i) {
      int b = (int) std::floor((buf[i] - p.gmin) / bw);
      if (b < 0) b = 0; if (b >= nb) b = nb - 1;
      h[b] += 1.0;
    }
    double e = 0.0;
    for (int b = 0; b < nb; ++b)
      if (h[b] > 0.0) { const double q = h[b] / k; e -= q * std::log(q); }
    ent = e;
    if (ss_needs(ids, 12, 12)) {
      double q = 0.0;
      for (int a = 0; a < nb; ++a) {
        if (h[a] == 0.0) continue;
        const double pa = h[a] / k, ca = p.gmin + (a + 0.5) * bw;
        for (int b = 0; b < nb; ++b) {
          if (h[b] == 0.0) continue;
          q += std::fabs(ca - (p.gmin + (b + 0.5) * bw)) * pa * (h[b] / k);
        }
      }
      rao = q;
    }
  }

  /* ---- lag-1 autocorrelation, slope, interfacial area --------------- */
  double r1 = NA_REAL, sdq = NA_REAL, sdr = NA_REAL;
  if (ss_needs(ids, 20, 22)) {
    double g = 0.0; int ng = 0;
    for (int r = 0; r < wr; ++r)
      for (int c = 0; c + 1 < wc; ++c) {
        const double d = buf[r * wc + c] - buf[r * wc + c + 1];
        g += d * d; ++ng;
      }
    for (int r = 0; r + 1 < wr; ++r)
      for (int c = 0; c < wc; ++c) {
        const double d = buf[r * wc + c] - buf[(r + 1) * wc + c];
        g += d * d; ++ng;
      }
    const double gamma1 = (ng > 0) ? g / (2.0 * ng) : NA_REAL;
    if (m2 > 0.0) r1 = 1.0 - gamma1 / m2;

    if (ss_needs(ids, 21, 22)) {
      double ssl = 0.0, sar = 0.0; int nd = 0;
      for (int r = 1; r + 1 < wr; ++r)
        for (int c = 1; c + 1 < wc; ++c) {
          const double gx = (buf[r * wc + c + 1] - buf[r * wc + c - 1]) /
                            (2.0 * p.resx);
          const double gy = (buf[(r - 1) * wc + c] - buf[(r + 1) * wc + c]) /
                            (2.0 * p.resy);
          const double s2 = gx * gx + gy * gy;
          ssl += s2; sar += std::sqrt(1.0 + s2) - 1.0; ++nd;
        }
      if (nd > 0) { sdq = std::sqrt(ssl / nd); sdr = sar / nd; }
    }
  }

  /* ---- grey-scale gliding-box lacunarity ---------------------------- */
  double lac = NA_REAL;
  if (ss_needs(ids, 23, 23)) {
    const int b = p.lacbox;
    if (b >= 1 && b <= wr && b <= wc) {
      double sm = 0.0, sm2 = 0.0; int nbx = 0;
      for (int r = 0; r + b <= wr; ++r)
        for (int c = 0; c + b <= wc; ++c) {
          double mass = 0.0;
          for (int rr = r; rr < r + b; ++rr)
            for (int cc = c; cc < c + b; ++cc)
              mass += buf[rr * wc + cc] - p.gmin;
          sm += mass; sm2 += mass * mass; ++nbx;
        }
      if (nbx > 0) {
        const double em = sm / nbx;
        if (em > 0.0) lac = (sm2 / nbx) / (em * em);
      }
    }
  }

  /* ---- grey level co-occurrence matrix ------------------------------ */
  double gcon = NA_REAL, gdis = NA_REAL, ghom = NA_REAL,
         gasm = NA_REAL, gent = NA_REAL, gcor = NA_REAL;
  if (need_glcm) {
    const int L = p.nlevels;
    const double bw = (p.gmax > p.gmin) ? (p.gmax - p.gmin) / L : 1.0;
    std::vector<int> q(k);
    for (int i = 0; i < k; ++i) {
      int b = (int) std::floor((buf[i] - p.gmin) / bw);
      if (b < 0) b = 0; if (b >= L) b = L - 1;
      q[i] = b;
    }
    std::vector<double> P((size_t) L * L, 0.0);
    double npair = 0.0;
    const int dr[4] = {0, -1, -1, -1};
    const int dc[4] = {1,  1,  0, -1};
    for (int d = 0; d < 4; ++d)
      for (int r = 0; r < wr; ++r)
        for (int c = 0; c < wc; ++c) {
          const int r2 = r + dr[d], c2 = c + dc[d];
          if (r2 < 0 || r2 >= wr || c2 < 0 || c2 >= wc) continue;
          const int a = q[r * wc + c], bq = q[r2 * wc + c2];
          P[(size_t) a * L + bq] += 1.0;      /* symmetric GLCM */
          P[(size_t) bq * L + a] += 1.0;
          npair += 2.0;
        }
    if (npair > 0.0) {
      for (size_t i = 0; i < P.size(); ++i) P[i] /= npair;
      double con = 0, dis = 0, hom = 0, asm_ = 0, ent2 = 0;
      double mu = 0, mu2 = 0, cov = 0;
      for (int a = 0; a < L; ++a)
        for (int b = 0; b < L; ++b) {
          const double pv = P[(size_t) a * L + b];
          if (pv <= 0.0) continue;
          const double diff = (double) (a - b);
          con  += pv * diff * diff;
          dis  += pv * std::fabs(diff);
          hom  += pv / (1.0 + diff * diff);
          asm_ += pv * pv;
          ent2 -= pv * std::log(pv);
          mu   += pv * a;
        }
      for (int a = 0; a < L; ++a)
        for (int b = 0; b < L; ++b) {
          const double pv = P[(size_t) a * L + b];
          if (pv <= 0.0) continue;
          mu2 += pv * (a - mu) * (a - mu);
          cov += pv * (a - mu) * (b - mu);
        }
      gcon = con; gdis = dis; ghom = hom; gasm = asm_; gent = ent2;
      gcor = (mu2 > 0.0) ? cov / mu2 : NA_REAL;
    }
  }

  for (size_t j = 0; j < ids.size(); ++j) {
    double v = NA_REAL;
    switch (ids[j]) {
      case 1:  v = mean; break;
      case 2:  v = sad;  break;
      case 3:  v = sq;   break;
      case 4:  v = (m2 > 0.0) ? m3 / (sq * sq * sq) : NA_REAL; break;
      case 5:  v = (m2 > 0.0) ? m4 / (m2 * m2) : NA_REAL; break;
      case 6:  v = mx - mean; break;
      case 7:  v = mean - mn; break;
      case 8:  v = mx - mn;   break;
      case 9:  v = s10z; break;
      case 10: v = (mean != 0.0) ? sq / std::fabs(mean) : NA_REAL; break;
      case 11: v = ent;  break;
      case 12: v = rao;  break;
      case 20: v = r1;   break;
      case 21: v = sdq;  break;
      case 22: v = sdr;  break;
      case 23: v = lac;  break;
      case 30: v = gcon; break;
      case 31: v = gdis; break;
      case 32: v = ghom; break;
      case 33: v = gasm; break;
      case 34: v = gent; break;
      case 35: v = gcor; break;
    }
    out[j] = v;
  }
}

/* rank transfer inside a window: keep the arrangement of `ranks_of`,
 * take the value multiset of `values_of`                                */
inline void ss_rank_transfer(const double* ranks_of, const double* values_of,
                             int k, std::vector<double>& sorted,
                             std::vector<int>& ord, double* out) {
  sorted.assign(values_of, values_of + k);
  std::sort(sorted.begin(), sorted.end());
  ord.resize(k);
  for (int i = 0; i < k; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) {
              if (ranks_of[a] == ranks_of[b]) return a < b;
              return ranks_of[a] < ranks_of[b];
            });
  for (int i = 0; i < k; ++i) out[ord[i]] = sorted[i];
}

#endif
