#pragma once
#include <string>
#include <vector>
#include <unordered_map>

// Struct to hold adjacent marker pairs and their LD values for efficient access in C++.
struct LDPair {
  std::string Name1;
  std::string Name2;
  double      LD;
};

// Main function to create blocks from LD data, called from R via Rcpp.
std::vector<std::vector<std::string>> make_blocks(
    const std::unordered_map<std::string, double>& ld_lookup,
    std::vector<LDPair>                           ld_adj,
    const std::vector<std::string>&                marker_names,
    const std::unordered_map<std::string, int>&    marker_idx,
    std::unordered_map<std::string, bool>&         assigned,
    const std::string&                             method,
    double                                         threshold,
    int                                            tolerance,
    bool                                           tol_reset,
    const std::string&                             start
);
