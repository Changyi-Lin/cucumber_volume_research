#pragma once

#include "types.hpp"

#include <array>
#include <filesystem>
#include <string>
#include <vector>

enum class DenoiseMethod { none, sor, ror, sor_ror };

struct DenoiseConfig {
  DenoiseMethod method = DenoiseMethod::none;
  int sor_k = 30;
  double sor_std_ratio = 2.0;
  double ror_radius_mm = 2.5;
  int ror_min_neighbors = 3;
  // 0 disables protection. Otherwise this is a fraction of PCA length per end.
  double end_protection_fraction = 0.0;
  // Region used for metrics even when protection is disabled.
  double metric_end_fraction = 0.05;
  bool clustering = false;
  double cluster_radius_mm = 3.0;
  bool component_keep_largest = true;
  std::size_t component_min_size = 1;
  // A detached end component is preserved only when it is spatially close to the main body.
  double component_end_keep_distance_multiplier = 2.0;
  // Set only after human/spatial inspection confirms that every removed
  // disconnected component is noise.  Its length contribution is then
  // reported separately instead of being mislabeled as end damage.
  bool confirmed_disconnected_noise = false;

  bool adaptive_local = false;
  int local_knn = 20;
  int local_pca_knn = 20;
  double local_mad_factor = 5.0;
  bool local_regional_threshold = true;
  double local_end_threshold_multiplier = 1.5;
  double local_density_weight = 0.5;
  double local_plane_weight = 0.3;
  double local_variation_weight = 0.2;
  double local_score_threshold = 0.8;

  // Expensive vector copies are enabled only for explicit stage inspection runs.
  bool capture_stages = false;
};

struct DistributionStats {
  double mean = 0, median = 0, std = 0, p5 = 0, p25 = 0, p75 = 0, p95 = 0;
};

struct DensityAnalysis {
  DistributionStats overall, bottom, middle, top;
  std::size_t overall_points = 0, bottom_points = 0, middle_points = 0, top_points = 0;
  double runtime_ms = 0;
};

struct RobustStats {
  double median = 0, mad = 0, threshold = 0;
};

struct LocalRegionStats {
  RobustStats density, plane_residual, surface_variation;
};

struct ComponentInfo {
  std::size_t id = 0, point_count = 0;
  std::array<double, 3> centroid{};
  std::array<double, 3> min{};
  std::array<double, 3> max{};
  double distance_to_main_mm = 0;
  bool main_component = false;
  bool kept = false;
};

struct DenoiseMetrics {
  std::size_t points_before = 0, points_after = 0, removed_points = 0;
  double removed_percent = 0;
  double top_removed_percent = 0, middle_removed_percent = 0, bottom_removed_percent = 0;
  std::size_t top_before = 0, middle_before = 0, bottom_before = 0;
  std::size_t top_after = 0, middle_after = 0, bottom_after = 0;
  double length_before_mm = 0, length_after_mm = 0, length_change_mm = 0;
  double length_change_percent = 0;
  double confirmed_component_length_change_mm = 0;
  double unexplained_length_change_mm = 0;
  bool possible_over_denoise = false;
  bool warning_end_damage = false;
  bool invalid_over_filtering = false;
  double sor_ms = 0, ror_ms = 0, clustering_ms = 0, local_ms = 0, denoise_total_ms = 0;
  std::size_t number_of_clusters = 0, largest_cluster_points = 0;
  std::size_t second_cluster_points = 0, removed_cluster_points = 0;
  std::size_t points_after_sor = 0, points_after_ror = 0;
  std::size_t points_after_component = 0, points_after_local = 0;
  std::size_t removed_sor_points = 0, removed_ror_points = 0;
  std::size_t removed_local_points = 0;
  LocalRegionStats local_bottom, local_middle, local_top;
};

struct DenoiseStages {
  std::vector<Point> after_sor, removed_sor;
  std::vector<Point> after_ror, removed_ror;
  std::vector<Point> after_component, removed_component;
  std::vector<Point> after_local, removed_local;
};

struct DenoiseOutput {
  std::vector<Point> points;
  std::vector<Point> removed;
  DenoiseMetrics metrics;
  std::vector<ComponentInfo> components;
  std::vector<std::size_t> component_labels;
  DenoiseStages stages;
};

DenoiseMethod parse_denoise_method(const std::string& name);
std::string denoise_method_name(DenoiseMethod method);

DenoiseOutput denoise_points(const std::vector<Point>& input,
                             const DenoiseConfig& config,
                             const std::vector<Point>* protection_reference = nullptr);

DensityAnalysis analyze_nearest_neighbor_density(const std::vector<Point>& points,
                                                 double end_fraction = 0.05);

void write_xyz_ply(const std::filesystem::path& path, const std::vector<Point>& points);
void write_component_rgb_ply(const std::filesystem::path& path,
                             const std::vector<Point>& points,
                             const std::vector<std::size_t>& labels,
                             const std::vector<ComponentInfo>& components);
void write_component_csv(const std::filesystem::path& path,
                         const std::vector<ComponentInfo>& components);
