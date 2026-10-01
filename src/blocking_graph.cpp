// Compute-bound routines for the graph-based haploblocking in R/blocking_graph.R.
//
// The LD calculation these used to sit beside now lives in src/pairwise_ld.cpp,
// since it is not specific to this method. What remains is: one LD question that
// is asked differently from a pairwise table - whether any member of a block
// reaches a threshold against one marker, which can stop at the first member
// that does - and the connected-component search that reads blocks off the graph.

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <functional>
using namespace Rcpp;


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
