#pragma once

#include "types.hpp"

#include <vector>

ComponentCleanupResult cleanup_tiny_components(Mesh& mesh,
                                               double minimum_largest_face_ratio = 0.95);

ValidationResult validate_mesh(const Mesh& mesh, const std::vector<Point>& input,
                               const DatasetStats& input_stats,
                               double offset_mm);
