// Compiled version of ld_func() in R/pairwise_LD.R. Computes the same r^2 as
// cor(use = "pairwise.complete.obs"), over the samples observed at both markers.

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

//' pairwise_ld_cpp
//'
//' Pairwise r^2 between markers on one chromosome. Called by ld_func_c().
//'
//' geno    : samples in rows, markers in columns
//' window  : compare markers at most this many positions apart; pass the marker
//'           count for every pair
//' min_r2  : drop pairs below this r^2; pass a negative value to keep all
//' min_obs : samples that must be observed at both markers for a pair to be kept
//'
//' Returns Locus1, Locus2 (1-based marker positions) and LD. Pairs with an
//' undefined r^2 are left out.
// [[Rcpp::export]]
DataFrame pairwise_ld_cpp(NumericMatrix geno,
                          int window,
                          double min_r2,
                          int min_obs) {
  int n_sample = geno.nrow();
  int n_snp = geno.ncol();

  std::vector<int> idxA;
  std::vector<int> idxB;
  std::vector<double> R2;

  // Reserve room for the most pairs the window allows, capped to limit memory
  double reserve_pairs = (double) n_snp * std::min(window, n_snp);
  int reserve_n = (int) std::min(reserve_pairs, 5e7);
  idxA.reserve(std::max(1000, reserve_n));
  idxB.reserve(std::max(1000, reserve_n));
  R2.reserve(std::max(1000, reserve_n));

  for (int i = 0; i < n_snp - 1; i++) {
    int j_end = std::min(n_snp - 1, i + window);

    for (int j = i + 1; j <= j_end; j++) {
      double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
      int n = 0;

      for (int k = 0; k < n_sample; k++) {
        double x = geno(k, i);
        double y = geno(k, j);
        if (!NumericVector::is_na(x) && !NumericVector::is_na(y)) {
          sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y; n++;
        }
      }

      double r2 = NA_REAL;

      // r^2 is undefined when either marker is constant across the shared samples
      if (n >= min_obs) {
        double num = n * sxy - sx * sy;
        double den_x = n * sxx - sx * sx;
        double den_y = n * syy - sy * sy;
        if (den_x > 0.0 && den_y > 0.0) {
          double r = num / std::sqrt(den_x * den_y);
          r2 = r * r;
        }
      }

      // Skip undefined pairs and those below min_r2
      if (std::isnan(r2) || r2 < min_r2) continue;

      idxA.push_back(i + 1);
      idxB.push_back(j + 1);
      R2.push_back(r2);
    }
  }

  return DataFrame::create(
    Named("Locus1") = idxA,
    Named("Locus2") = idxB,
    Named("LD") = R2
  );
}
