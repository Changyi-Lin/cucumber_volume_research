#pragma once

#include "types.hpp"

#include <filesystem>
#include <vector>

std::vector<Point> read_ply_xyz(const std::filesystem::path& path);
DatasetStats compute_dataset_stats(const std::vector<Point>& points);

