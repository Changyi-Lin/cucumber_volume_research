#include "voxel_downsample.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>

namespace {
struct Accumulator {
  double x = 0, y = 0, z = 0;
  std::uint32_t count = 0;
};
}

std::vector<Point> voxel_downsample_centroid(const std::vector<Point>& points,
                                             const double voxel_size_mm,
                                             const DatasetStats& stats) {
  if (!(voxel_size_mm > 0)) return points;
  constexpr std::uint64_t mask = (std::uint64_t{1} << 21) - 1;
  const auto dim = [&](int k) {
    return static_cast<std::uint64_t>(std::floor((stats.max[k] - stats.min[k]) / voxel_size_mm)) + 1;
  };
  if (dim(0) > mask || dim(1) > mask || dim(2) > mask)
    throw std::runtime_error("Voxel grid exceeds 21-bit packed-key range");

  std::unordered_map<std::uint64_t, Accumulator> voxels;
  voxels.reserve(std::min<std::size_t>(points.size() / 3 + 1, 800000));
  const double inv = 1.0 / voxel_size_mm;
  for (const auto& p : points) {
    const auto ix = static_cast<std::uint64_t>(std::floor((p.x() - stats.min[0]) * inv));
    const auto iy = static_cast<std::uint64_t>(std::floor((p.y() - stats.min[1]) * inv));
    const auto iz = static_cast<std::uint64_t>(std::floor((p.z() - stats.min[2]) * inv));
    const std::uint64_t key = ix | (iy << 21) | (iz << 42);
    auto& a = voxels[key];
    a.x += p.x(); a.y += p.y(); a.z += p.z(); ++a.count;
  }
  std::vector<Point> out;
  out.reserve(voxels.size());
  for (const auto& [key, a] : voxels) {
    (void)key;
    const double inv_count = 1.0 / static_cast<double>(a.count);
    out.emplace_back(a.x * inv_count, a.y * inv_count, a.z * inv_count);
  }
  return out;
}

std::vector<Point> systematic_sample(const std::vector<Point>& points,
                                     const std::size_t target_count) {
  if (target_count == 0 || target_count >= points.size()) return points;
  std::vector<Point> out;
  out.reserve(target_count);
  const long double step = static_cast<long double>(points.size()) / target_count;
  for (std::size_t i = 0; i < target_count; ++i) {
    const auto index = std::min<std::size_t>(static_cast<std::size_t>((i + 0.5L) * step), points.size() - 1);
    out.push_back(points[index]);
  }
  return out;
}

