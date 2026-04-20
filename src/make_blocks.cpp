#include <Rcpp.h>
#include <algorithm>
#include <unordered_map>
#include "extend_block.h"
#include "make_blocks.h"
using namespace Rcpp;


// Drives block formation for a single chromosome. Finds seed marker pairs and
// calls extend_block to grow each block left and right. Mirrors the R
// implementation — see make_blocks in R/def_haploblocks.R.
//
// ld_lookup   : "Name1,Name2" -> LD value
// ld_adj      : adjacent marker pairs with Name1, Name2, LD fields (ld_adj in R)
// marker_names: ordered list of marker names for the chromosome
// marker_idx  : marker name -> 0-based position index
// assigned    : marker name -> already in a block (mutated in place)
// method      : "flanking" or "average"
// threshold   : minimum LD value to seed or extend a block
// tolerance   : consecutive below-threshold markers allowed during extension
// tol_reset   : if true, reset tolerance counter after each accepted marker
// start       : "LD" (highest-LD pair first) or "beginning" (left to right sweep)
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
) {

  std::vector<std::vector<std::string>> chrom_blocks;
  int n_markers = (int)marker_names.size();

  if(start == "LD") {

    while(true) {

      auto   max_it    = std::max_element(ld_adj.begin(), ld_adj.end(), [](const LDPair& a, const LDPair& b) { return a.LD < b.LD; });
      int    max_ld_idx = (int)(max_it - ld_adj.begin());
      double max_ld     = ld_adj[max_ld_idx].LD;
      if(max_ld <= threshold) break;

      std::string seed_marker_left  = ld_adj[max_ld_idx].Name1;
      std::string seed_marker_right = ld_adj[max_ld_idx].Name2;

      // if either seed marker is already in a block, drop this pair and move on
      if(assigned[seed_marker_left] || assigned[seed_marker_right]) {
        ld_adj.erase(ld_adj.begin() + max_ld_idx);
        continue;
      }

      std::vector<std::string> block = {seed_marker_left, seed_marker_right};

      // Perform block extensions in both directions from the seed pair
      block = extend_block(-1, seed_marker_left,  marker_names, marker_idx, assigned, ld_lookup, block, method, threshold, tolerance, tol_reset);
      block = extend_block( 1, seed_marker_right, marker_names, marker_idx, assigned, ld_lookup, block, method, threshold, tolerance, tol_reset);

      // Mark each block as assigned
      for(const std::string& marker : block) {
        assigned[marker] = true;
      }
      // Add to final output list of blocks for this chromosome
      chrom_blocks.push_back(block);

      // Returns true if either Name1 or Name2 of the LDPair is already assigned to a block
      auto pair_in_block = [&](const LDPair& p) {
        return assigned[p.Name1] || assigned[p.Name2];
      };

      // Erase all adjacent pairs from ld_adj that involve any marker already in a block
      ld_adj.erase(std::remove_if(ld_adj.begin(), ld_adj.end(), pair_in_block), ld_adj.end());

      // If no adjacent pairs remain, break the loop
      if(ld_adj.empty()) break;
    }
  
  } else if(start == "beginning") {  // beginning

    int position   = 0;
    int total_snps = n_markers;

    // Loop through markers from left to right, using each unassigned marker as a seed to grow a block
    while(true) {
      
      // The first unassigned marker becomes the seed for a new block
      const std::string& seed_marker = marker_names[position];
      std::vector<std::string> block = {seed_marker};

      block = extend_block(1, seed_marker, marker_names, marker_idx, assigned, ld_lookup, block, method, threshold, tolerance, tol_reset);

      // Mark each block as assigned
      for(const std::string& marker : block) {
        assigned[marker] = true;
      }
      // Add to final output list of blocks for this chromosome
      chrom_blocks.push_back(block);

      // Move position to the next unassigned marker after the current block
      position = marker_idx.at(block.back()) + 1;
      // If we've reached the end of the marker list, break the loop
      if(position >= total_snps) break;
    }
  }

  return chrom_blocks;
}


//' make_blocks_cpp
//'
//' R-facing wrapper around make_blocks. Converts R named vectors to C++ types
//' once up front, calls make_blocks, then converts the result back to R types.
//' See make_blocks in R/def_haploblocks.R for full parameter documentation.
// [[Rcpp::export]]
List make_blocks_cpp(
    NumericVector   ld_lookup_r,
    DataFrame       ld_adj_r,
    CharacterVector marker_names_r,
    IntegerVector   marker_idx_r,
    LogicalVector   assigned_r,
    std::string     method,
    double          threshold,
    int             tolerance,
    bool            tol_reset,
    std::string     start
) {

  // Build ld_lookup map from R named numeric vector: "Name1,Name2" -> LD value
  std::unordered_map<std::string, double> ld_lookup;
  {
    CharacterVector names = ld_lookup_r.names();
    int n = ld_lookup_r.size();
    ld_lookup.reserve(n);

    // Build the ld_lookup map from the R named numeric vector. 
    for(int i = 0; i < n; i++) {
      ld_lookup[as<std::string>(names[i])] = ld_lookup_r[i];
    }
  }

  // Build 0-based — R's marker_idx uses seq_along (1-based) but C++ vectors are 0-based
  std::unordered_map<std::string, int> marker_idx;
  {
    CharacterVector names = marker_idx_r.names();
    int n = marker_idx_r.size();
    marker_idx.reserve(n);

    // Build the marker_idx map from the R named integer vector, converting to 0-based indexing for C++.
    for(int i = 0; i < n; i++) {
      marker_idx[as<std::string>(names[i])] = marker_idx_r[i] - 1;
    }
  }

  std::unordered_map<std::string, bool> assigned;
  {
    CharacterVector names = assigned_r.names();
    int n = assigned_r.size();
    assigned.reserve(n);

    // Build the assigned map from the R named logical vector.
    for(int i = 0; i < n; i++) {
      assigned[as<std::string>(names[i])] = (bool)assigned_r[i];
    }
  }

  // Convert marker_names from R CharacterVector to C++ vector<string>
  std::vector<std::string> marker_names = std::vector<std::string>(marker_names_r.begin(), marker_names_r.end());

  // Convert ld_adj from R DataFrame to C++ vector<LDPair>
  CharacterVector col_name1 = ld_adj_r["Name1"];
  CharacterVector col_name2 = ld_adj_r["Name2"];
  NumericVector   col_ld    = ld_adj_r["LD"];
  int n_pairs = col_ld.size();

  // Reserve space to avoid multiple reallocations during push_back
  std::vector<LDPair> ld_pairs;
  ld_pairs.reserve(n_pairs);

  // Populate ld_pairs with LDPair structs built from the R DataFrame columns
  for(int i = 0; i < n_pairs; i++) {
    ld_pairs.push_back({as<std::string>(col_name1[i]), as<std::string>(col_name2[i]), col_ld[i]});
  }

  // Call make_blocks with the converted C++ types to get the list of blocks for this chromosome
  std::vector<std::vector<std::string>> chrom_blocks = make_blocks(
    ld_lookup, ld_pairs, marker_names, marker_idx, assigned,
    method, threshold, tolerance, tol_reset, start
  );

  // Convert chrom_blocks back to an R list of character vectors
  List blocks_out(chrom_blocks.size());
  for(int i = 0; i < (int)chrom_blocks.size(); i++) {
    // Wrap converts vector<string> to CharacterVector, which can be stored in an R list
    blocks_out[i] = wrap(chrom_blocks[i]);
  }

  // Rebuild assigned output vector from the updated map
  // Create memory for the output LogicalVector
  int n_markers = (int)marker_names.size();
  LogicalVector assigned_out(n_markers);

  // Set names for the output LogicalVector to match the input marker names
  assigned_out.names() = marker_names_r;
  // Populate the output LogicalVector with values from the assigned map, ensuring the order matches marker_names
  for(int i = 0; i < n_markers; i++) {
    assigned_out[i] = assigned[marker_names[i]];
  }

  // Return the final list of blocks
  return List::create(blocks_out, assigned_out);
}
