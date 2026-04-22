#include "extend_block.h"
#include <cmath>
#include <numeric>

// Internal C++ extend_block called directly from make_blocks_cpp without
// crossing the R boundary on each iteration.
std::vector<std::string> extend_block(
    int direction,
    std::string edge_marker,
    const std::vector<std::string>& marker_names,
    const std::unordered_map<std::string, int>& marker_idx,
    std::unordered_map<std::string, bool>& assigned,
    const std::unordered_map<std::string, double>& ld_lookup,
    std::vector<std::string> block,
    const std::string& method,
    double threshold,
    int tolerance,
    bool tol_reset
) {

  // Counter for consecutive failures to meet LD threshold
  int tolerance_counter = 0;

  // List of failed candidate markers not yet added to a block
  std::vector<std::string> failed_markers;

  // The total number of markers in the dataset (length of marker_names)
  int n_markers = (int)marker_names.size();

  // Loop until we reach a break condition
  while(true) {
    int edge_index = marker_idx.at(edge_marker);

    // Check if we've reached the end of the marker list in the given direction
    if(direction == -1 && edge_index == 0)             break;
    if(direction ==  1 && edge_index == n_markers - 1) break;

    int step = (int)failed_markers.size() + 1;
    // Calculate the candidate marker index based on the current edge and direction
    int candidate_index = (direction == -1) ? (edge_index - step) : (edge_index + step);

    // If the candidate index is out of bounds, break the loop
    if(candidate_index < 0 || candidate_index >= n_markers) break;

    // Get the candidate marker name from the marker_names vector
    const std::string& candidate_marker = marker_names[candidate_index];

    // If the candidate marker is already assigned to a block, break the loop
    if(assigned.count(candidate_marker) && assigned.at(candidate_marker)) break;

    // Vector to hold LD values for the candidate marker
    std::vector<double> ld_vals;

    // Flanking method: LD is calculated between the candidate marker and the edge marker only
    if(method == "flanking") {
      if(direction == -1) {
        auto it = ld_lookup.find(candidate_marker + "," + edge_marker);
        // Only insert if the LD value exists and is not NaN
        if(it != ld_lookup.end() && !std::isnan(it->second)) {
          ld_vals.push_back(it->second);
        }
      } else {
        // candidate is to the right so Name1 = edge, Name2 = candidate
        auto it = ld_lookup.find(edge_marker + "," + candidate_marker);
        // Only insert if the LD value exists and is not NaN
        if(it != ld_lookup.end() && !std::isnan(it->second)) {
          ld_vals.push_back(it->second);
        }
      }
    }

    // Average method: LD is calculated between the candidate marker and all markers in the block, then averaged
    if(method == "average") {
      if(direction == -1) {
        // candidate is to the left so Name1 = candidate, Name2 = any block marker
        for(const auto& bm : block) {
          auto it = ld_lookup.find(candidate_marker + "," + bm);
          // Only insert if the LD value exists and is not NaN
          if(it != ld_lookup.end() && !std::isnan(it->second)) {
            ld_vals.push_back(it->second);
          }
        }
      } else {
        // candidate is to the right so Name1 = any block marker, Name2 = candidate
        for(const auto& bm : block) {
          auto it = ld_lookup.find(bm + "," + candidate_marker);
          // Only insert if the LD value exists and is not NaN
          if(it != ld_lookup.end() && !std::isnan(it->second)) {
            ld_vals.push_back(it->second);
          }
        }
      }
    }

    // No LD data found for this candidate — can't evaluate it, treat as below threshold
    if(ld_vals.empty()) {
      if(tolerance_counter >= tolerance) break;
      tolerance_counter++;
      failed_markers.insert(failed_markers.begin(), candidate_marker);
      continue;
    }

    // Calculate the mean LD value for the candidate marker, C++ has no built-in mean function so we sum and divide by count
    double ld_mean = std::accumulate(ld_vals.begin(), ld_vals.end(), 0.0) / ld_vals.size();

    // Check if the mean LD meets the threshold
    if(ld_mean <= threshold) {
      // If it doesn't meet the threshold, check if we've reached the tolerance limit
      if(tolerance_counter >= tolerance) {
        break;
      } else {
        // If we haven't reached the tolerance limit, increment the counter and add the candidate marker to the list of failed markers
        tolerance_counter++;
        failed_markers.insert(failed_markers.begin(), candidate_marker);
        continue;
      }
    }

    // Reset the tolerance counter if tol_reset is true
    if(tol_reset) tolerance_counter = 0;

    
    if(direction == -1) {
      // Prepend to block: [candidate | failed_markers | block]
      block.insert(block.begin(), failed_markers.begin(), failed_markers.end());
      block.insert(block.begin(), candidate_marker);
    } else {
      // Append to block: [block | failed_markers | candidate]
      block.insert(block.end(), failed_markers.begin(), failed_markers.end());
      block.push_back(candidate_marker);
    }

    // Clear the list of failed markers now
    failed_markers.clear();

    // Edge marker is now the candidate marker for the next iteration
    edge_marker = candidate_marker;
  }

  return block;
}
