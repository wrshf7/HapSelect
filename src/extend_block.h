#pragma once
#include <string>
#include <vector>
#include <unordered_map>

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
);
