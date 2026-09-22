// Compute-bound routines for the graph-based haploblocking in R/blocking_graph.R.
//
// compute_local_ld_edges_cpp and has_strong_ld_to_block_cpp both calculate
// pairwise-complete Pearson r^2 over a sample-by-marker genotype matrix. They are
// separate because they are asked different questions: the first builds the sparse
// LD edge table for a whole chromosome within a forward marker window, while the
// second answers one yes/no question about a single marker against one block, and
// so can stop at the first member that clears the threshold.
//
// Both skip a pair with fewer than three samples observed at both markers, and a
// pair where either marker has no variance among those samples. The three-sample
// floor matters: r is exactly +/-1 for two complete observations, which would
// otherwise plant perfect LD on a pair with no evidence behind it.

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <functional>
using namespace Rcpp;


//' compute_local_ld_edges_cpp
//'
//' Windowed pairwise r^2 over one chromosome. geno is sample by marker; snps and
//' pos follow the same column order. For marker i only i+1, ..., i+W_snp are
//' evaluated, and only pairs reaching min_r2 are returned, with indices converted
//' back to 1-based R positions. See local_ld_edges in R/blocking_graph.R.
// [[Rcpp::export]]
DataFrame compute_local_ld_edges_cpp(NumericMatrix geno,
                                     CharacterVector snps,
                                     NumericVector pos,
                                     int W_snp,
                                     double min_r2) {
  int n_sample = geno.nrow();
  int n_snp = geno.ncol();

  std::vector<std::string> SNP_A;
  std::vector<std::string> SNP_B;
  std::vector<int> idxA;
  std::vector<int> idxB;
  std::vector<double> bpA;
  std::vector<double> bpB;
  std::vector<double> R2;

  int reserve_n = std::max(1000, n_snp * 5);
  SNP_A.reserve(reserve_n);
  SNP_B.reserve(reserve_n);
  idxA.reserve(reserve_n);
  idxB.reserve(reserve_n);
  bpA.reserve(reserve_n);
  bpB.reserve(reserve_n);
  R2.reserve(reserve_n);

  for (int i = 0; i < n_snp - 1; i++) {
    int j_end = std::min(n_snp - 1, i + W_snp);

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

      if (n < 3) continue;
      double num = n * sxy - sx * sy;
      double den_x = n * sxx - sx * sx;
      double den_y = n * syy - sy * sy;
      if (den_x <= 0.0 || den_y <= 0.0) continue;

      double r = num / std::sqrt(den_x * den_y);
      double r2 = r * r;
      if (!std::isnan(r2) && r2 >= min_r2) {
        SNP_A.push_back(as<std::string>(snps[i]));
        SNP_B.push_back(as<std::string>(snps[j]));
        idxA.push_back(i + 1);
        idxB.push_back(j + 1);
        bpA.push_back(pos[i]);
        bpB.push_back(pos[j]);
        R2.push_back(r2);
      }
    }
  }

  return DataFrame::create(
    Named("SNP_A") = SNP_A, Named("SNP_B") = SNP_B,
    Named("idxA") = idxA, Named("idxB") = idxB,
    Named("bpA") = bpA, Named("bpB") = bpB, Named("R2") = R2
  );
}

//' has_strong_ld_to_block_cpp
//'
//' Whether a target marker reaches threshold r^2 against any one member of a
//' block, recomputed from the genotypes rather than read from the edge table, so
//' neither the LD window nor the edge-table floor constrains it. Returns as soon
//' as a member qualifies. See block_ld_support in R/blocking_graph.R.
// [[Rcpp::export]]
bool has_strong_ld_to_block_cpp(NumericMatrix geno,
                                int target_col,
                                IntegerVector member_cols,
                                double threshold) {
  int n_sample = geno.nrow();
  int t = target_col - 1;
  double threshold2 = threshold;

  for (int m = 0; m < member_cols.size(); m++) {
    int j = member_cols[m] - 1;
    if (j < 0 || j >= geno.ncol() || j == t) continue;

    double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
    int n = 0;

    for (int k = 0; k < n_sample; k++) {
      double x = geno(k, t);
      double y = geno(k, j);
      if (!NumericVector::is_na(x) && !NumericVector::is_na(y)) {
        sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y; n++;
      }
    }

    if (n < 3) continue;
    double num = n * sxy - sx * sy;
    double den_x = n * sxx - sx * sx;
    double den_y = n * syy - sy * sy;
    if (den_x <= 0.0 || den_y <= 0.0) continue;

    double r = num / std::sqrt(den_x * den_y);
    double r2 = r * r;
    if (!std::isnan(r2) && r2 >= threshold2) return true;
  }

  return false;
}

//' connected_components_cpp
//'
//' Component label for every vertex of an undirected graph, by union-find with
//' path compression and union by size. Stands in for igraph's components() so
//' that the one graph operation the method needs does not pull in igraph.
//'
//' Labels are assigned in order of first appearance while walking the vertices,
//' which is not meaningful on its own: core_blocks() re-sorts each component by
//' marker index and re-orders the components by their first marker, so the
//' labelling here cannot affect the blocks that come out.
//'
//' from, to : 1-based vertex indices, one pair per edge
//' n_vertex : total vertices, including any with no edges at all
// [[Rcpp::export]]
IntegerVector connected_components_cpp(IntegerVector from, IntegerVector to,
                                       int n_vertex) {
  std::vector<int> parent(n_vertex);
  std::vector<int> size(n_vertex, 1);
  for (int i = 0; i < n_vertex; i++) parent[i] = i;

  // Find with full path compression
  std::function<int(int)> root = [&](int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };

  for (int e = 0; e < from.size(); e++) {
    int a = root(from[e] - 1);
    int b = root(to[e] - 1);
    if (a == b) continue;
    // Attach the smaller tree to the larger, to keep the trees shallow
    if (size[a] < size[b]) std::swap(a, b);
    parent[b] = a;
    size[a] += size[b];
  }

  // Relabel roots to 1, 2, ... in order of first appearance
  IntegerVector membership(n_vertex);
  std::vector<int> label(n_vertex, 0);
  int next_label = 0;
  for (int v = 0; v < n_vertex; v++) {
    int r = root(v);
    if (label[r] == 0) label[r] = ++next_label;
    membership[v] = label[r];
  }

  return membership;
}
