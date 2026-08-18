#pragma once

#include "types.hpp"

#include <vector>

std::vector<Point> voxel_downsample_centroid(const std::vector<Point>& points,
                                             double voxel_size_mm,
                                             const DatasetStats& stats);
std::vector<Point> systematic_sample(const std::vector<Point>& points,
                                     std::size_t target_count);

