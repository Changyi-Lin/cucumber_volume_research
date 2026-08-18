#include "denoise.hpp"

#include "ply_reader.hpp"

#include <CGAL/Kd_tree.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace {

using Clock = std::chrono::steady_clock;
using SearchTraits = CGAL::Search_traits_3<Kernel>;
using NeighborSearch = CGAL::Orthogonal_k_neighbor_search<SearchTraits>;
using SearchTree = NeighborSearch::Tree;

double elapsed_ms(const Clock::time_point begin, const Clock::time_point end) {
  return std::chrono::duration<double, std::milli>(end - begin).count();
}

struct ProjectionReference {
  std::array<double, 3> axis{};
  double min = 0, max = 0;
};

double projection(const Point& p, const ProjectionReference& ref) {
  return p.x() * ref.axis[0] + p.y() * ref.axis[1] + p.z() * ref.axis[2];
}

ProjectionReference make_projection_reference(const std::vector<Point>& points) {
  const DatasetStats stats = compute_dataset_stats(points);
  ProjectionReference ref;
  ref.axis = stats.pca_axes[2];
  ref.min = std::numeric_limits<double>::infinity();
  ref.max = -std::numeric_limits<double>::infinity();
  for (const auto& point : points) {
    const double s = projection(point, ref);
    ref.min = std::min(ref.min, s);
    ref.max = std::max(ref.max, s);
  }
  return ref;
}

bool protected_end(const Point& p, const ProjectionReference& ref, const double fraction) {
  if (!(fraction > 0)) return false;
  const double s = projection(p, ref);
  const double margin = fraction * (ref.max - ref.min);
  return s <= ref.min + margin || s >= ref.max - margin;
}

std::pair<std::vector<Point>, std::vector<Point>> sor_filter(
    const std::vector<Point>& input, const int k, const double std_ratio,
    const ProjectionReference& ref, const double protection_fraction, double& runtime_ms) {
  const auto begin = Clock::now();
  if (input.size() <= static_cast<std::size_t>(k + 1))
    throw std::runtime_error("SOR input has fewer points than K+1");
  SearchTree tree(input.begin(), input.end());
  tree.build();
  std::vector<double> mean_distances(input.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (long long i = 0; i < static_cast<long long>(input.size()); ++i) {
    NeighborSearch search(tree, input[static_cast<std::size_t>(i)], k + 1);
    double sum = 0;
    int count = 0;
    for (auto it = search.begin(); it != search.end(); ++it) {
      const double distance = std::sqrt(std::max(0.0, CGAL::to_double(it->second)));
      if (distance <= 1e-12) continue;
      sum += distance;
      if (++count == k) break;
    }
    mean_distances[static_cast<std::size_t>(i)] = count ? sum / count : 0;
  }
  const double mean = std::accumulate(mean_distances.begin(), mean_distances.end(), 0.0) /
                      static_cast<double>(mean_distances.size());
  double variance = 0;
  for (const double value : mean_distances) variance += (value - mean) * (value - mean);
  variance /= static_cast<double>(mean_distances.size());
  const double threshold = mean + std_ratio * std::sqrt(variance);
  std::vector<Point> kept, removed;
  kept.reserve(input.size());
  removed.reserve(input.size() / 50 + 1);
  for (std::size_t i = 0; i < input.size(); ++i) {
    if (mean_distances[i] <= threshold || protected_end(input[i], ref, protection_fraction))
      kept.push_back(input[i]);
    else
      removed.push_back(input[i]);
  }
  runtime_ms = elapsed_ms(begin, Clock::now());
  return {std::move(kept), std::move(removed)};
}

std::pair<std::vector<Point>, std::vector<Point>> ror_filter(
    const std::vector<Point>& input, const double radius, const int min_neighbors,
    const ProjectionReference& ref, const double protection_fraction, double& runtime_ms) {
  const auto begin = Clock::now();
  if (input.size() <= static_cast<std::size_t>(min_neighbors + 1))
    throw std::runtime_error("ROR input has fewer points than min_neighbors+1");
  SearchTree tree(input.begin(), input.end());
  tree.build();
  std::vector<unsigned char> keep(input.size(), 0);
  const double radius_sq = radius * radius;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (long long i = 0; i < static_cast<long long>(input.size()); ++i) {
    const auto index = static_cast<std::size_t>(i);
    if (protected_end(input[index], ref, protection_fraction)) {
      keep[index] = 1;
      continue;
    }
    NeighborSearch search(tree, input[index], min_neighbors + 1);
    int within = 0;
    for (auto it = search.begin(); it != search.end(); ++it) {
      const double distance_sq = CGAL::to_double(it->second);
      if (distance_sq <= 1e-24) continue;
      if (distance_sq <= radius_sq) ++within;
    }
    keep[index] = within >= min_neighbors;
  }
  std::vector<Point> kept, removed;
  kept.reserve(input.size());
  removed.reserve(input.size() / 50 + 1);
  for (std::size_t i = 0; i < input.size(); ++i) {
    (keep[i] ? kept : removed).push_back(input[i]);
  }
  runtime_ms = elapsed_ms(begin, Clock::now());
  return {std::move(kept), std::move(removed)};
}

struct GridKeyHash {
  std::size_t operator()(const std::array<int, 3>& key) const noexcept {
    std::uint64_t h = 1469598103934665603ULL;
    for (int v : key) { h ^= static_cast<std::uint32_t>(v); h *= 1099511628211ULL; }
    return static_cast<std::size_t>(h);
  }
};

struct ClusterOutput {
  std::vector<std::size_t> labels;
  std::vector<std::size_t> sizes;
};

ClusterOutput radius_clusters(const std::vector<Point>& points, const double radius) {
  std::unordered_map<std::array<int, 3>, std::vector<std::size_t>, GridKeyHash> grid;
  grid.reserve(points.size() / 2 + 1);
  const double inv = 1.0 / radius;
  const auto key_of = [&](const Point& p) {
    return std::array<int, 3>{{static_cast<int>(std::floor(p.x() * inv)),
                               static_cast<int>(std::floor(p.y() * inv)),
                               static_cast<int>(std::floor(p.z() * inv))}};
  };
  for (std::size_t i = 0; i < points.size(); ++i) grid[key_of(points[i])].push_back(i);
  const double radius_sq = radius * radius;
  ClusterOutput output;
  output.labels.assign(points.size(), std::numeric_limits<std::size_t>::max());
  std::queue<std::size_t> queue;
  for (std::size_t start = 0; start < points.size(); ++start) {
    if (output.labels[start] != std::numeric_limits<std::size_t>::max()) continue;
    const std::size_t label = output.sizes.size();
    output.sizes.push_back(0);
    output.labels[start] = label;
    queue.push(start);
    while (!queue.empty()) {
      const std::size_t current = queue.front(); queue.pop();
      ++output.sizes[label];
      const auto base = key_of(points[current]);
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
          for (int dz = -1; dz <= 1; ++dz) {
            const std::array<int, 3> key{{base[0] + dx, base[1] + dy, base[2] + dz}};
            const auto bucket = grid.find(key);
            if (bucket == grid.end()) continue;
            for (const std::size_t neighbor : bucket->second) {
              if (output.labels[neighbor] != std::numeric_limits<std::size_t>::max()) continue;
              const double x = points[current].x() - points[neighbor].x();
              const double y = points[current].y() - points[neighbor].y();
              const double z = points[current].z() - points[neighbor].z();
              if (x * x + y * y + z * z <= radius_sq) {
                output.labels[neighbor] = label;
                queue.push(neighbor);
              }
            }
          }
    }
  }
  return output;
}

DistributionStats distribution(std::vector<double> values) {
  DistributionStats result;
  if (values.empty()) return result;
  std::sort(values.begin(), values.end());
  const auto quantile = [&](double q) {
    const double index = q * (values.size() - 1);
    const auto lo = static_cast<std::size_t>(std::floor(index));
    const auto hi = static_cast<std::size_t>(std::ceil(index));
    const double t = index - lo;
    return values[lo] * (1 - t) + values[hi] * t;
  };
  result.mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
  double variance = 0;
  for (const double value : values) variance += (value - result.mean) * (value - result.mean);
  result.std = std::sqrt(variance / values.size());
  result.p5 = quantile(0.05); result.p25 = quantile(0.25);
  result.median = quantile(0.50); result.p75 = quantile(0.75); result.p95 = quantile(0.95);
  return result;
}

RobustStats robust_stats(std::vector<double> values, const double factor) {
  RobustStats result;
  if (values.empty()) return result;
  std::sort(values.begin(), values.end());
  const auto median_of_sorted = [](const std::vector<double>& sorted) {
    const std::size_t middle = sorted.size() / 2;
    return sorted.size() % 2 ? sorted[middle] : 0.5 * (sorted[middle - 1] + sorted[middle]);
  };
  result.median = median_of_sorted(values);
  std::vector<double> deviations;
  deviations.reserve(values.size());
  for (const double value : values) deviations.push_back(std::abs(value - result.median));
  std::sort(deviations.begin(), deviations.end());
  result.mad = median_of_sorted(deviations);
  result.threshold = result.median + factor * result.mad;
  return result;
}

int metric_region(const Point& point, const ProjectionReference& ref, const double fraction) {
  const double margin = fraction * (ref.max - ref.min);
  const double s = projection(point, ref);
  if (s <= ref.min + margin) return 0;
  if (s >= ref.max - margin) return 2;
  return 1;
}

struct LocalFeatures {
  double density = 0;
  double plane_residual = 0;
  double surface_variation = 0;
};

std::vector<LocalFeatures> compute_local_features(const std::vector<Point>& points,
                                                  const int density_k,
                                                  const int pca_k) {
  const int query_k = std::max(density_k, pca_k);
  if (points.size() <= static_cast<std::size_t>(query_k + 1))
    throw std::runtime_error("Adaptive local input has fewer points than K+1");
  SearchTree tree(points.begin(), points.end());
  tree.build();
  std::vector<LocalFeatures> features(points.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (long long raw_index = 0; raw_index < static_cast<long long>(points.size()); ++raw_index) {
    const std::size_t index = static_cast<std::size_t>(raw_index);
    NeighborSearch search(tree, points[index], query_k + 1);
    std::vector<Eigen::Vector3d> neighbors;
    neighbors.reserve(static_cast<std::size_t>(query_k));
    double density_sum = 0;
    int density_count = 0;
    for (auto it = search.begin(); it != search.end(); ++it) {
      const double distance_sq = std::max(0.0, CGAL::to_double(it->second));
      if (distance_sq <= 1e-24) continue;
      const Point neighbor = it->first;
      neighbors.emplace_back(neighbor.x(), neighbor.y(), neighbor.z());
      if (density_count < density_k) {
        density_sum += std::sqrt(distance_sq);
        ++density_count;
      }
      if (static_cast<int>(neighbors.size()) == query_k) break;
    }
    LocalFeatures feature;
    feature.density = density_count ? density_sum / density_count : 0;
    const int used = std::min<int>(pca_k, static_cast<int>(neighbors.size()));
    if (used >= 3) {
      Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
      for (int i = 0; i < used; ++i) centroid += neighbors[static_cast<std::size_t>(i)];
      centroid /= used;
      Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
      for (int i = 0; i < used; ++i) {
        const Eigen::Vector3d delta = neighbors[static_cast<std::size_t>(i)] - centroid;
        covariance.noalias() += delta * delta.transpose();
      }
      covariance /= used;
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance, Eigen::ComputeEigenvectors);
      if (solver.info() == Eigen::Success) {
        const Eigen::Vector3d eigenvalues = solver.eigenvalues().cwiseMax(0.0);
        const Eigen::Vector3d normal = solver.eigenvectors().col(0);
        const Eigen::Vector3d center(points[index].x(), points[index].y(), points[index].z());
        feature.plane_residual = std::abs(normal.dot(center - centroid));
        const double sum = eigenvalues.sum();
        feature.surface_variation = sum > 0 ? eigenvalues[0] / sum : 0;
      }
    }
    features[index] = feature;
  }
  return features;
}

struct LocalFilterOutput {
  std::vector<Point> kept, removed;
  LocalRegionStats bottom, middle, top;
};

LocalFilterOutput adaptive_local_filter(const std::vector<Point>& points,
                                        const DenoiseConfig& config,
                                        const ProjectionReference& ref) {
  const std::vector<LocalFeatures> features = compute_local_features(
      points, config.local_knn, config.local_pca_knn);
  const double region_fraction = config.end_protection_fraction > 0
      ? config.end_protection_fraction : config.metric_end_fraction;
  std::array<std::vector<double>, 3> density, plane, variation;
  for (std::size_t i = 0; i < points.size(); ++i) {
    const int region = metric_region(points[i], ref, region_fraction);
    density[region].push_back(features[i].density);
    plane[region].push_back(features[i].plane_residual);
    variation[region].push_back(features[i].surface_variation);
  }
  if (!config.local_regional_threshold) {
    std::vector<double> all_density, all_plane, all_variation;
    for (int region = 0; region < 3; ++region) {
      all_density.insert(all_density.end(), density[region].begin(), density[region].end());
      all_plane.insert(all_plane.end(), plane[region].begin(), plane[region].end());
      all_variation.insert(all_variation.end(), variation[region].begin(), variation[region].end());
    }
    for (int region = 0; region < 3; ++region) {
      density[region] = all_density;
      plane[region] = all_plane;
      variation[region] = all_variation;
    }
  }
  std::array<LocalRegionStats, 3> region_stats;
  for (int region = 0; region < 3; ++region) {
    const double factor = config.local_mad_factor *
        (region == 1 ? 1.0 : config.local_end_threshold_multiplier);
    region_stats[region].density = robust_stats(density[region], factor);
    region_stats[region].plane_residual = robust_stats(plane[region], factor);
    region_stats[region].surface_variation = robust_stats(variation[region], factor);
  }
  const auto normalized = [](const double value, const RobustStats& stats) {
    const double scale = std::max(stats.threshold - stats.median,
                                  std::max(1e-12, std::abs(stats.median) * 1e-9));
    return std::max(0.0, (value - stats.median) / scale);
  };
  LocalFilterOutput output;
  output.kept.reserve(points.size());
  output.removed.reserve(points.size() / 100 + 1);
  for (std::size_t i = 0; i < points.size(); ++i) {
    const int region = metric_region(points[i], ref, region_fraction);
    const auto& stats = region_stats[region];
    const double density_score = normalized(features[i].density, stats.density);
    const double plane_score = normalized(features[i].plane_residual, stats.plane_residual);
    const double variation_score = normalized(features[i].surface_variation, stats.surface_variation);
    const double score = config.local_density_weight * density_score +
                         config.local_plane_weight * plane_score +
                         config.local_variation_weight * variation_score;
    const bool low_density = density_score > 1.0;
    const bool inconsistent_surface = plane_score > 1.0 || variation_score > 1.0;
    const bool remove = low_density && inconsistent_surface &&
                        score >= config.local_score_threshold;
    (remove ? output.removed : output.kept).push_back(points[i]);
  }
  output.bottom = region_stats[0];
  output.middle = region_stats[1];
  output.top = region_stats[2];
  return output;
}

std::vector<ComponentInfo> describe_components(const std::vector<Point>& points,
                                                const ClusterOutput& clusters,
                                                const std::vector<std::size_t>& order) {
  std::vector<ComponentInfo> descriptions(clusters.sizes.size());
  for (std::size_t label = 0; label < descriptions.size(); ++label) {
    auto& info = descriptions[label];
    info.id = label;
    info.point_count = clusters.sizes[label];
    info.min.fill(std::numeric_limits<double>::infinity());
    info.max.fill(-std::numeric_limits<double>::infinity());
  }
  for (std::size_t i = 0; i < points.size(); ++i) {
    auto& info = descriptions[clusters.labels[i]];
    const std::array<double, 3> coordinates{{points[i].x(), points[i].y(), points[i].z()}};
    for (int axis = 0; axis < 3; ++axis) {
      info.centroid[axis] += coordinates[axis];
      info.min[axis] = std::min(info.min[axis], coordinates[axis]);
      info.max[axis] = std::max(info.max[axis], coordinates[axis]);
    }
  }
  for (auto& info : descriptions)
    for (double& coordinate : info.centroid) coordinate /= info.point_count;
  if (order.empty()) return descriptions;
  const std::size_t main_label = order[0];
  descriptions[main_label].main_component = true;
  std::vector<Point> main_points;
  main_points.reserve(clusters.sizes[main_label]);
  for (std::size_t i = 0; i < points.size(); ++i)
    if (clusters.labels[i] == main_label) main_points.push_back(points[i]);
  SearchTree main_tree(main_points.begin(), main_points.end());
  main_tree.build();
  for (std::size_t label = 0; label < descriptions.size(); ++label) {
    if (label == main_label) continue;
    double minimum_sq = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < points.size(); ++i) {
      if (clusters.labels[i] != label) continue;
      NeighborSearch search(main_tree, points[i], 1);
      if (search.begin() != search.end())
        minimum_sq = std::min(minimum_sq, CGAL::to_double(search.begin()->second));
    }
    descriptions[label].distance_to_main_mm = std::sqrt(std::max(0.0, minimum_sq));
  }
  return descriptions;
}

}  // namespace

DenoiseMethod parse_denoise_method(const std::string& input_name) {
  std::string name = input_name;
  std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) { return std::tolower(c); });
  if (name == "none" || name == "voxel") return DenoiseMethod::none;
  if (name == "sor") return DenoiseMethod::sor;
  if (name == "ror") return DenoiseMethod::ror;
  if (name == "sor_ror" || name == "sor+ror") return DenoiseMethod::sor_ror;
  throw std::runtime_error("Unknown denoise method: " + name);
}

std::string denoise_method_name(const DenoiseMethod method) {
  switch (method) {
    case DenoiseMethod::none: return "NONE";
    case DenoiseMethod::sor: return "SOR";
    case DenoiseMethod::ror: return "ROR";
    case DenoiseMethod::sor_ror: return "SOR_ROR";
  }
  return "UNKNOWN";
}

DenoiseOutput denoise_points(const std::vector<Point>& input, const DenoiseConfig& config,
                             const std::vector<Point>* protection_reference) {
  const auto total_begin = Clock::now();
  DenoiseOutput output;
  output.metrics.points_before = input.size();
  const auto& reference_points = protection_reference ? *protection_reference : input;
  const ProjectionReference ref = make_projection_reference(reference_points);
  std::vector<Point> current = input;

  if (config.method == DenoiseMethod::sor || config.method == DenoiseMethod::sor_ror) {
    auto [kept, removed] = sor_filter(current, config.sor_k, config.sor_std_ratio, ref,
                                      config.end_protection_fraction, output.metrics.sor_ms);
    output.metrics.removed_sor_points = removed.size();
    if (config.capture_stages) output.stages.removed_sor = removed;
    output.removed.insert(output.removed.end(), removed.begin(), removed.end());
    current = std::move(kept);
  }
  output.metrics.points_after_sor = current.size();
  if (config.capture_stages) output.stages.after_sor = current;

  if (config.method == DenoiseMethod::ror || config.method == DenoiseMethod::sor_ror) {
    auto [kept, removed] = ror_filter(current, config.ror_radius_mm, config.ror_min_neighbors, ref,
                                      config.end_protection_fraction, output.metrics.ror_ms);
    output.metrics.removed_ror_points = removed.size();
    if (config.capture_stages) output.stages.removed_ror = removed;
    output.removed.insert(output.removed.end(), removed.begin(), removed.end());
    current = std::move(kept);
  }
  output.metrics.points_after_ror = current.size();
  if (config.capture_stages) output.stages.after_ror = current;

  const auto projected_length = [&](const std::vector<Point>& points) {
    if (points.empty()) return 0.0;
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    for (const auto& point : points) {
      const double value = projection(point, ref);
      minimum = std::min(minimum, value);
      maximum = std::max(maximum, value);
    }
    return maximum - minimum;
  };
  const double length_before_component = projected_length(current);

  if (config.clustering && !current.empty()) {
    const auto cluster_begin = Clock::now();
    ClusterOutput clusters = radius_clusters(current, config.cluster_radius_mm);
    std::vector<std::size_t> order(clusters.sizes.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
      return clusters.sizes[a] > clusters.sizes[b];
    });
    output.metrics.number_of_clusters = order.size();
    if (!order.empty()) output.metrics.largest_cluster_points = clusters.sizes[order[0]];
    if (order.size() > 1) output.metrics.second_cluster_points = clusters.sizes[order[1]];
    output.components = describe_components(current, clusters, order);
    output.component_labels = clusters.labels;
    const bool safe_to_remove = order.size() > 1 &&
        static_cast<double>(output.metrics.largest_cluster_points) / current.size() >= 0.98 &&
        static_cast<double>(output.metrics.second_cluster_points) / current.size() < 0.01;
    if (safe_to_remove) {
      std::vector<Point> kept;
      kept.reserve(output.metrics.largest_cluster_points);
      std::vector<Point> removed;
      removed.reserve(current.size() - output.metrics.largest_cluster_points);
      const std::size_t main_label = order[0];
      std::vector<unsigned char> keep_component(order.size(), 0);
      keep_component[main_label] = 1;
      for (std::size_t label = 0; label < output.components.size(); ++label) {
        if (label == main_label) continue;
        bool touches_end = false;
        for (std::size_t i = 0; i < current.size(); ++i) {
          if (clusters.labels[i] == label &&
              protected_end(current[i], ref, config.end_protection_fraction)) {
            touches_end = true;
            break;
          }
        }
        const bool plausible_detached_end = touches_end &&
            output.components[label].distance_to_main_mm <=
                config.component_end_keep_distance_multiplier * config.cluster_radius_mm &&
            output.components[label].point_count >= config.component_min_size;
        const bool meets_minimum = !config.component_keep_largest &&
            output.components[label].point_count >= config.component_min_size;
        keep_component[label] = plausible_detached_end || meets_minimum;
      }
      for (std::size_t i = 0; i < current.size(); ++i) {
        if (keep_component[clusters.labels[i]])
          kept.push_back(current[i]);
        else
          removed.push_back(current[i]);
      }
      for (std::size_t label = 0; label < output.components.size(); ++label)
        output.components[label].kept = keep_component[label];
      output.metrics.removed_cluster_points = removed.size();
      if (config.capture_stages) output.stages.removed_component = removed;
      output.removed.insert(output.removed.end(), removed.begin(), removed.end());
      current = std::move(kept);
    } else {
      for (auto& component : output.components) component.kept = true;
    }
    output.metrics.clustering_ms = elapsed_ms(cluster_begin, Clock::now());
  }
  output.metrics.points_after_component = current.size();
  if (config.confirmed_disconnected_noise) {
    output.metrics.confirmed_component_length_change_mm = std::max(
        0.0, length_before_component - projected_length(current));
  }
  if (config.capture_stages) output.stages.after_component = current;

  if (config.adaptive_local && !current.empty()) {
    const auto local_begin = Clock::now();
    LocalFilterOutput local = adaptive_local_filter(current, config, ref);
    output.metrics.removed_local_points = local.removed.size();
    output.metrics.local_bottom = local.bottom;
    output.metrics.local_middle = local.middle;
    output.metrics.local_top = local.top;
    if (config.capture_stages) output.stages.removed_local = local.removed;
    output.removed.insert(output.removed.end(), local.removed.begin(), local.removed.end());
    current = std::move(local.kept);
    output.metrics.local_ms = elapsed_ms(local_begin, Clock::now());
  }
  output.metrics.points_after_local = current.size();
  if (config.capture_stages) output.stages.after_local = current;

  output.points = std::move(current);
  output.metrics.points_after = output.points.size();
  output.metrics.removed_points = input.size() - output.points.size();
  output.metrics.removed_percent = 100.0 * output.metrics.removed_points / input.size();

  const auto region = [&](const Point& point) {
    return metric_region(point, ref, config.metric_end_fraction);
  };
  for (const auto& point : input) {
    const int r = region(point);
    if (r == 0) ++output.metrics.bottom_before;
    else if (r == 2) ++output.metrics.top_before;
    else ++output.metrics.middle_before;
  }
  for (const auto& point : output.points) {
    const int r = region(point);
    if (r == 0) ++output.metrics.bottom_after;
    else if (r == 2) ++output.metrics.top_after;
    else ++output.metrics.middle_after;
  }
  const auto removed_percent = [](std::size_t before, std::size_t after) {
    return before ? 100.0 * static_cast<double>(before - after) / before : 0.0;
  };
  output.metrics.bottom_removed_percent = removed_percent(output.metrics.bottom_before, output.metrics.bottom_after);
  output.metrics.middle_removed_percent = removed_percent(output.metrics.middle_before, output.metrics.middle_after);
  output.metrics.top_removed_percent = removed_percent(output.metrics.top_before, output.metrics.top_after);
  output.metrics.length_before_mm = ref.max - ref.min;
  double after_min = std::numeric_limits<double>::infinity();
  double after_max = -std::numeric_limits<double>::infinity();
  for (const auto& point : output.points) {
    const double s = projection(point, ref);
    after_min = std::min(after_min, s); after_max = std::max(after_max, s);
  }
  output.metrics.length_after_mm = output.points.empty() ? 0 : after_max - after_min;
  output.metrics.length_change_mm = output.metrics.length_before_mm - output.metrics.length_after_mm;
  output.metrics.length_change_percent = 100.0 * output.metrics.length_change_mm /
                                         output.metrics.length_before_mm;
  output.metrics.unexplained_length_change_mm = std::max(
      0.0, output.metrics.length_change_mm -
               output.metrics.confirmed_component_length_change_mm);
  output.metrics.warning_end_damage = output.metrics.unexplained_length_change_mm > 1.0 ||
      100.0 * output.metrics.unexplained_length_change_mm /
          output.metrics.length_before_mm > 0.5;
  output.metrics.invalid_over_filtering =
      output.metrics.unexplained_length_change_mm > 2.0;
  output.metrics.possible_over_denoise = output.metrics.warning_end_damage;
  output.metrics.denoise_total_ms = elapsed_ms(total_begin, Clock::now());
  return output;
}

DensityAnalysis analyze_nearest_neighbor_density(const std::vector<Point>& points,
                                                 const double end_fraction) {
  const auto begin = Clock::now();
  DensityAnalysis result;
  if (points.size() < 2) return result;
  const ProjectionReference ref = make_projection_reference(points);
  SearchTree tree(points.begin(), points.end()); tree.build();
  std::vector<double> all(points.size(), 0), bottom, middle, top;
  const double margin = end_fraction * (ref.max - ref.min);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (long long i = 0; i < static_cast<long long>(points.size()); ++i) {
    NeighborSearch search(tree, points[static_cast<std::size_t>(i)], 2);
    double nearest = 0;
    for (auto it = search.begin(); it != search.end(); ++it) {
      const double distance = std::sqrt(std::max(0.0, CGAL::to_double(it->second)));
      if (distance > 1e-12) { nearest = distance; break; }
    }
    all[static_cast<std::size_t>(i)] = nearest;
  }
  for (std::size_t i = 0; i < points.size(); ++i) {
    const double s = projection(points[i], ref);
    if (s <= ref.min + margin) bottom.push_back(all[i]);
    else if (s >= ref.max - margin) top.push_back(all[i]);
    else middle.push_back(all[i]);
  }
  result.overall = distribution(all); result.bottom = distribution(bottom);
  result.middle = distribution(middle); result.top = distribution(top);
  result.overall_points = all.size(); result.bottom_points = bottom.size();
  result.middle_points = middle.size(); result.top_points = top.size();
  result.runtime_ms = elapsed_ms(begin, Clock::now());
  return result;
}

void write_xyz_ply(const std::filesystem::path& path, const std::vector<Point>& points) {
  if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("Cannot write point PLY: " + path.string());
  out << "ply\nformat binary_little_endian 1.0\ncomment cucumber denoising benchmark\n"
      << "element vertex " << points.size() << "\nproperty double x\nproperty double y\n"
      << "property double z\nend_header\n";
  for (const auto& p : points) {
    const double xyz[3] = {p.x(), p.y(), p.z()};
    out.write(reinterpret_cast<const char*>(xyz), sizeof(xyz));
  }
}

void write_component_rgb_ply(const std::filesystem::path& path,
                             const std::vector<Point>& points,
                             const std::vector<std::size_t>& labels,
                             const std::vector<ComponentInfo>& components) {
  if (points.size() != labels.size())
    throw std::runtime_error("Component PLY point/label size mismatch");
  if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("Cannot write component PLY: " + path.string());
  out << "ply\nformat binary_little_endian 1.0\ncomment component labels: main=white removed=colored\n"
      << "element vertex " << points.size() << "\nproperty double x\nproperty double y\n"
      << "property double z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nend_header\n";
  static constexpr std::array<std::array<unsigned char, 3>, 8> palette{{
      {{230, 57, 70}}, {{46, 196, 182}}, {{69, 123, 157}}, {{255, 183, 3}},
      {{131, 56, 236}}, {{251, 133, 0}}, {{0, 150, 199}}, {{214, 40, 40}}
  }};
  for (std::size_t i = 0; i < points.size(); ++i) {
    const Point& point = points[i];
    const double xyz[3] = {point.x(), point.y(), point.z()};
    out.write(reinterpret_cast<const char*>(xyz), sizeof(xyz));
    const std::size_t label = labels[i];
    std::array<unsigned char, 3> color{{220, 220, 220}};
    if (label < components.size() && !components[label].main_component) {
      color = components[label].kept ? std::array<unsigned char, 3>{{160, 160, 160}}
                                     : palette[label % palette.size()];
    }
    out.write(reinterpret_cast<const char*>(color.data()), 3);
  }
}

void write_component_csv(const std::filesystem::path& path,
                         const std::vector<ComponentInfo>& components) {
  if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Cannot write component CSV: " + path.string());
  out << "component_id,point_count,centroid_x,centroid_y,centroid_z,"
         "min_x,min_y,min_z,max_x,max_y,max_z,distance_to_main_mm,main_component,kept\n";
  for (const auto& component : components) {
    out << std::setprecision(12) << component.id << ',' << component.point_count << ','
        << component.centroid[0] << ',' << component.centroid[1] << ',' << component.centroid[2] << ','
        << component.min[0] << ',' << component.min[1] << ',' << component.min[2] << ','
        << component.max[0] << ',' << component.max[1] << ',' << component.max[2] << ','
        << component.distance_to_main_mm << ',' << component.main_component << ','
        << component.kept << '\n';
  }
}
