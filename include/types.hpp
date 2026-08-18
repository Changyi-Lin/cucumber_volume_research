#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <array>
#include <cstddef>
#include <limits>
#include <string>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;

struct DatasetStats {
  std::size_t point_count = 0;
  std::array<double, 3> min{{0, 0, 0}};
  std::array<double, 3> max{{0, 0, 0}};
  std::array<double, 3> centroid{{0, 0, 0}};
  // Columns, ordered from smallest to largest covariance eigenvalue.
  std::array<std::array<double, 3>, 3> pca_axes{};
  std::array<double, 3> eigenvalues{{0, 0, 0}};
  std::array<double, 3> pca_extents{{0, 0, 0}};
  double aabb_diagonal = 0;
};

struct ValidationResult {
  std::size_t vertices = 0;
  std::size_t faces = 0;
  std::size_t edges = 0;
  std::size_t connected_components = 0;
  bool is_closed = false;
  bool is_triangle_mesh = false;
  bool self_intersection = false;
  bool manifold = false;
  bool outward_oriented = false;
  bool two_sided_wrap = false;
  long long euler_characteristic = 0;
  double genus = std::numeric_limits<double>::quiet_NaN();
  double surface_area_mm2 = std::numeric_limits<double>::quiet_NaN();
  double median_line_intersections = std::numeric_limits<double>::quiet_NaN();
  std::size_t lines_with_four_or_more = 0;
  std::size_t tested_lines = 0;
};

struct ComponentCleanupResult {
  std::size_t raw_connected_components = 0;
  std::size_t removed_components = 0;
  double largest_component_face_ratio = 1.0;
};

struct BenchmarkResult {
  double voxel_size_mm = 0;
  std::size_t original_points = 0;
  std::size_t input_points = 0;
  double alpha_mm = 0;
  double offset_mm = 0;
  int run_index = 0;
  double ply_load_ms = 0;
  double downsample_ms = 0;
  double denoise_ms = 0;
  double alpha_wrap_ms = 0;
  double component_cleanup_ms = 0;
  double validation_ms = 0;
  double volume_ms = 0;
  double total_ms = 0;
  double reconstruction_only_ms = 0;
  double online_core_ms = 0;
  ComponentCleanupResult cleanup;
  ValidationResult validation;
  double signed_volume_mm3 = std::numeric_limits<double>::quiet_NaN();
  double volume_mm3 = std::numeric_limits<double>::quiet_NaN();
  bool valid = false;
  std::string failure_reason;
};
