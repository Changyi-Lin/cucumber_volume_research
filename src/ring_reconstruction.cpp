#include "ring_reconstruction.hpp"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;
using Vec3 = Eigen::Vector3d;
constexpr double kPi = 3.141592653589793238462643383279502884;

double elapsed_ms(const Clock::time_point begin, const Clock::time_point end) {
  return std::chrono::duration<double, std::milli>(end - begin).count();
}

Vec3 as_vec(const Point& p) { return Vec3(p.x(), p.y(), p.z()); }
Point as_point(const Vec3& p) { return Point(p.x(), p.y(), p.z()); }

struct Slice {
  double s = 0;
  std::vector<std::size_t> point_indices;
  Vec3 raw_center = Vec3::Zero();
  Vec3 center = Vec3::Zero();
  Vec3 tangent = Vec3::UnitZ();
  Vec3 u = Vec3::UnitX();
  Vec3 v = Vec3::UnitY();
  std::vector<double> radii;
  std::vector<double> original_radii;
  std::size_t populated_sectors = 0;
  double representative_radius = 0;
};

struct PcaFrame {
  Vec3 centroid = Vec3::Zero();
  Vec3 axis = Vec3::UnitZ();
  Vec3 transverse = Vec3::UnitX();
};

PcaFrame compute_pca(const std::vector<Point>& points) {
  PcaFrame result;
  for (const auto& p : points) result.centroid += as_vec(p);
  result.centroid /= static_cast<double>(points.size());
  Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
  for (const auto& p : points) {
    const Vec3 d = as_vec(p) - result.centroid;
    covariance.noalias() += d * d.transpose();
  }
  covariance /= static_cast<double>(points.size());
  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
  if (solver.info() != Eigen::Success) throw std::runtime_error("PCA eigensolver failed");
  result.transverse = solver.eigenvectors().col(0).normalized();
  result.axis = solver.eigenvectors().col(2).normalized();
  // Make runs deterministic despite the arbitrary eigenvector sign.
  Eigen::Index dominant = 0;
  result.axis.cwiseAbs().maxCoeff(&dominant);
  if (result.axis[dominant] < 0) result.axis = -result.axis;
  return result;
}

double quantile_sorted(const std::vector<double>& sorted, const double q) {
  if (sorted.empty()) return std::numeric_limits<double>::quiet_NaN();
  const double index = q * static_cast<double>(sorted.size() - 1);
  const auto lo = static_cast<std::size_t>(std::floor(index));
  const auto hi = static_cast<std::size_t>(std::ceil(index));
  const double fraction = index - static_cast<double>(lo);
  return sorted[lo] * (1.0 - fraction) + sorted[hi] * fraction;
}

void interpolate_empty_sectors(std::vector<double>& radii) {
  const std::size_t n = radii.size();
  std::vector<std::size_t> valid;
  for (std::size_t i = 0; i < n; ++i)
    if (std::isfinite(radii[i]) && radii[i] > 0) valid.push_back(i);
  if (valid.size() < 3) throw std::runtime_error("Ring has fewer than three populated angular sectors");
  for (std::size_t i = 0; i < n; ++i) {
    if (std::isfinite(radii[i]) && radii[i] > 0) continue;
    std::size_t backward = 1;
    while (backward < n && !std::isfinite(radii[(i + n - backward) % n])) ++backward;
    std::size_t forward = 1;
    while (forward < n && !std::isfinite(radii[(i + forward) % n])) ++forward;
    const double a = radii[(i + n - backward) % n];
    const double b = radii[(i + forward) % n];
    radii[i] = (a * static_cast<double>(forward) + b * static_cast<double>(backward)) /
               static_cast<double>(backward + forward);
  }
}

double median_copy(std::vector<double> values) {
  if (values.empty()) return 0;
  std::sort(values.begin(), values.end());
  return quantile_sorted(values, 0.5);
}

std::vector<double> filter_sector_radial_outliers(const std::vector<double>& sorted,
                                                   const double multiplier) {
  if (!(multiplier > 0) || sorted.size() < 4) return sorted;
  const double median = quantile_sorted(sorted, 0.5);
  std::vector<double> deviations;
  deviations.reserve(sorted.size());
  for (const double value : sorted) deviations.push_back(std::abs(value - median));
  std::sort(deviations.begin(), deviations.end());
  const double mad = quantile_sorted(deviations, 0.5);
  const double robust_sigma = 1.4826 * mad;
  if (robust_sigma < 1e-9) return sorted;
  const double threshold = multiplier * robust_sigma;
  std::vector<double> filtered;
  filtered.reserve(sorted.size());
  for (const double value : sorted)
    if (std::abs(value - median) <= threshold) filtered.push_back(value);
  // A sparse sector is more trustworthy unfiltered than accidentally emptied.
  return filtered.size() >= 3 ? filtered : sorted;
}

std::vector<double> smooth_ring_radius_circular(const std::vector<double>& radii,
                                                const std::size_t window) {
  if (window == 1) return radii;
  if ((window != 3 && window != 5 && window != 7) || window > radii.size())
    throw std::runtime_error("ring_smoothing_window must be 1, 3, 5, or 7");
  std::vector<double> output(radii.size());
  const long long half = static_cast<long long>(window / 2);
  for (std::size_t j = 0; j < radii.size(); ++j) {
    std::vector<double> neighborhood;
    neighborhood.reserve(window);
    for (long long d = -half; d <= half; ++d) {
      const auto index = static_cast<std::size_t>(
          (static_cast<long long>(j) + d + static_cast<long long>(radii.size())) %
          static_cast<long long>(radii.size()));
      neighborhood.push_back(radii[index]);
    }
    output[j] = median_copy(std::move(neighborhood));
  }
  return output;
}

RingReliabilityMetrics compute_ring_reliability(const std::vector<double>& radii,
                                                const std::size_t ring_index,
                                                const std::size_t point_count,
                                                const double coverage) {
  RingReliabilityMetrics metrics;
  metrics.ring_index = ring_index;
  metrics.point_count = point_count;
  metrics.angular_coverage = coverage;
  metrics.mean_radius = std::accumulate(radii.begin(), radii.end(), 0.0) /
                        static_cast<double>(radii.size());
  double variance = 0;
  double roughness = 0;
  for (std::size_t j = 0; j < radii.size(); ++j) {
    const double delta = radii[j] - metrics.mean_radius;
    variance += delta * delta;
    roughness += std::abs(radii[j] - radii[(j + radii.size() - 1) % radii.size()]);
  }
  metrics.radius_stddev = std::sqrt(variance / static_cast<double>(radii.size()));
  metrics.radius_cv = metrics.mean_radius > 1e-9
      ? metrics.radius_stddev / metrics.mean_radius : 0;
  metrics.angular_roughness = roughness / static_cast<double>(radii.size());
  metrics.normalized_roughness = metrics.mean_radius > 1e-9
      ? metrics.angular_roughness / metrics.mean_radius : 0;
  return metrics;
}

std::pair<double, double> robust_location_sigma(const std::vector<double>& values) {
  const double median = median_copy(values);
  std::vector<double> deviations;
  deviations.reserve(values.size());
  for (const double value : values) deviations.push_back(std::abs(value - median));
  return {median, 1.4826 * median_copy(std::move(deviations))};
}

struct EndpointSelection {
  std::vector<std::size_t> candidates;
  std::vector<std::size_t> shape_indices;
  std::vector<std::size_t> radius_indices;
  EndpointShapeTemplate shape_template;
  std::size_t boundary_index = 0;
};

EndpointSelection select_endpoint_reliable_rings(
    const std::vector<Slice>& rings, std::vector<RingReliabilityMetrics>& metrics,
    const bool bottom, const RingReconstructionConfig& config) {
  EndpointSelection selection;
  const std::size_t count = std::min(config.endpoint_candidate_rings, rings.size());
  for (std::size_t i = 0; i < count; ++i)
    selection.candidates.push_back(bottom ? i : rings.size() - 1 - i);
  std::vector<double> cvs, roughness, coverage, point_counts;
  for (const auto index : selection.candidates) {
    cvs.push_back(metrics[index].radius_cv);
    roughness.push_back(metrics[index].normalized_roughness);
    coverage.push_back(metrics[index].angular_coverage);
    point_counts.push_back(static_cast<double>(metrics[index].point_count));
  }
  const auto [cv_median, cv_sigma] = robust_location_sigma(cvs);
  const auto [rough_median, rough_sigma] = robust_location_sigma(roughness);
  const auto [coverage_median, coverage_sigma] = robust_location_sigma(coverage);
  const auto [points_median, points_sigma] = robust_location_sigma(point_counts);
  const double cv_threshold = std::min(config.endpoint_absolute_max_cv,
      cv_median + config.endpoint_quality_mad_multiplier * std::max(cv_sigma, 1e-6));
  const double rough_threshold = std::min(config.endpoint_absolute_max_normalized_roughness,
      rough_median + config.endpoint_quality_mad_multiplier * std::max(rough_sigma, 1e-6));
  const double coverage_threshold = std::max(0.5,
      coverage_median - config.endpoint_quality_mad_multiplier * coverage_sigma);
  const double point_threshold = std::max(12.0,
      points_median - config.endpoint_quality_mad_multiplier * points_sigma);
  for (const auto index : selection.candidates) {
    auto& metric = metrics[index];
    metric.reliable_for_shape = metric.angular_coverage >= coverage_threshold &&
        static_cast<double>(metric.point_count) >= point_threshold &&
        metric.radius_cv <= cv_threshold &&
        metric.normalized_roughness <= rough_threshold;
    metric.reliable_for_mesh_boundary = metric.reliable_for_shape;
    // Radius fitting is intentionally more permissive than shape fitting.
    metric.reliable_for_radius_fit = metric.angular_coverage >= 0.5 &&
        metric.point_count >= 12 && std::isfinite(metric.mean_radius) &&
        metric.mean_radius > 0;
    if (metric.reliable_for_shape &&
        selection.shape_indices.size() < config.endpoint_shape_ring_count)
      selection.shape_indices.push_back(index);
    if (metric.reliable_for_radius_fit &&
        selection.radius_indices.size() < config.endpoint_fit_rings)
      selection.radius_indices.push_back(index);
  }
  if (selection.shape_indices.empty()) {
    // Safe fallback: choose the least rough candidate instead of the outermost ring.
    const auto best = *std::min_element(selection.candidates.begin(), selection.candidates.end(),
        [&](const std::size_t a, const std::size_t b) {
          return metrics[a].radius_cv + metrics[a].normalized_roughness <
                 metrics[b].radius_cv + metrics[b].normalized_roughness;
        });
    selection.shape_indices.push_back(best);
    metrics[best].reliable_for_shape = true;
    metrics[best].reliable_for_mesh_boundary = true;
  }
  selection.boundary_index = selection.shape_indices.front();
  if (selection.radius_indices.size() < 4) {
    selection.radius_indices.clear();
    for (std::size_t i = 0; i < std::min<std::size_t>(config.endpoint_fit_rings,
                                                     selection.candidates.size()); ++i)
      selection.radius_indices.push_back(selection.candidates[i]);
  }
  selection.shape_template.source_ring_indices = selection.shape_indices;
  selection.shape_template.normalized_radii.resize(config.angular_bins);
  for (std::size_t sector = 0; sector < config.angular_bins; ++sector) {
    std::vector<double> normalized;
    for (const auto index : selection.shape_indices)
      normalized.push_back(rings[index].radii[sector] /
                           std::max(1e-9, metrics[index].mean_radius));
    selection.shape_template.normalized_radii[sector] = median_copy(std::move(normalized));
  }
  const double template_mean = std::accumulate(selection.shape_template.normalized_radii.begin(),
      selection.shape_template.normalized_radii.end(), 0.0) /
      static_cast<double>(selection.shape_template.normalized_radii.size());
  for (double& radius : selection.shape_template.normalized_radii)
    radius /= std::max(1e-9, template_mean);
  return selection;
}

std::vector<double> open_uniform_knots(const int control_count, const int degree) {
  std::vector<double> knots(static_cast<std::size_t>(control_count + degree + 1));
  for (int i = 0; i <= degree; ++i) knots[static_cast<std::size_t>(i)] = 0;
  const int interior = control_count - degree - 1;
  for (int i = 1; i <= interior; ++i)
    knots[static_cast<std::size_t>(degree + i)] =
        static_cast<double>(i) / static_cast<double>(interior + 1);
  for (int i = control_count; i < control_count + degree + 1; ++i)
    knots[static_cast<std::size_t>(i)] = 1;
  return knots;
}

double bspline_basis(const int i, const int degree, const double t,
                     const std::vector<double>& knots, const int control_count) {
  if (degree == 0) {
    if ((knots[static_cast<std::size_t>(i)] <= t &&
         t < knots[static_cast<std::size_t>(i + 1)]) ||
        (t == 1.0 && i == control_count - 1)) return 1.0;
    return 0.0;
  }
  double value = 0;
  const double left_den = knots[static_cast<std::size_t>(i + degree)] -
                          knots[static_cast<std::size_t>(i)];
  if (left_den > 0)
    value += (t - knots[static_cast<std::size_t>(i)]) / left_den *
             bspline_basis(i, degree - 1, t, knots, control_count);
  const double right_den = knots[static_cast<std::size_t>(i + degree + 1)] -
                           knots[static_cast<std::size_t>(i + 1)];
  if (right_den > 0)
    value += (knots[static_cast<std::size_t>(i + degree + 1)] - t) / right_den *
             bspline_basis(i + 1, degree - 1, t, knots, control_count);
  return value;
}

struct SplineProfile {
  int degree = 3;
  int control_count = 4;
  std::vector<double> knots;
  Eigen::VectorXd controls;
  double max_distance = 1;
  double boundary_distance = 0;
  double boundary_radius = 0;
};

std::vector<double> isotonic_non_decreasing(const std::vector<double>& input) {
  struct Block { double sum = 0; std::size_t count = 0; };
  std::vector<Block> blocks;
  for (const double value : input) {
    blocks.push_back({value, 1});
    while (blocks.size() >= 2) {
      const double previous = blocks[blocks.size() - 2].sum /
                              static_cast<double>(blocks[blocks.size() - 2].count);
      const double current = blocks.back().sum / static_cast<double>(blocks.back().count);
      if (previous <= current) break;
      blocks[blocks.size() - 2].sum += blocks.back().sum;
      blocks[blocks.size() - 2].count += blocks.back().count;
      blocks.pop_back();
    }
  }
  std::vector<double> output;
  output.reserve(input.size());
  for (const auto& block : blocks) {
    const double value = block.sum / static_cast<double>(block.count);
    for (std::size_t i = 0; i < block.count; ++i) output.push_back(value);
  }
  return output;
}

SplineProfile fit_endpoint(const std::vector<Slice>& rings, const bool bottom,
                           const std::vector<std::size_t>& ring_indices,
                           const double tip_s, const bool enforce_monotonic,
                           std::vector<RingReliabilityMetrics>* metrics) {
  const std::size_t k = ring_indices.size();
  if (k < 4) throw std::runtime_error("Cubic endpoint fitting requires at least four rings");
  std::vector<double> distance;
  std::vector<double> radius;
  distance.reserve(k + 1);
  radius.reserve(k + 1);
  distance.push_back(0);
  radius.push_back(0);
  std::vector<double> endpoint_radius;
  endpoint_radius.reserve(k);
  for (const auto index : ring_indices) endpoint_radius.push_back(rings[index].representative_radius);
  for (std::size_t i = 1; i < endpoint_radius.size(); ++i) {
    if (endpoint_radius[i] + 1e-6 < endpoint_radius[i - 1] && metrics) {
      (*metrics)[ring_indices[i - 1]].endpoint_radius_rebound = true;
      (*metrics)[ring_indices[i]].endpoint_radius_rebound = true;
    }
  }
  if (enforce_monotonic) endpoint_radius = isotonic_non_decreasing(endpoint_radius);
  for (std::size_t i = 0; i < k; ++i) {
    const Slice& ring = rings[ring_indices[i]];
    distance.push_back(std::abs(ring.s - tip_s));
    radius.push_back(endpoint_radius[i]);
  }
  std::vector<std::size_t> order(distance.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](const std::size_t a, const std::size_t b) {
    return distance[a] < distance[b];
  });
  SplineProfile result;
  result.max_distance = std::max(1e-9, *std::max_element(distance.begin(), distance.end()));
  const Slice& boundary = bottom ? rings.front() : rings.back();
  result.boundary_distance = std::abs(boundary.s - tip_s);
  result.boundary_radius = boundary.representative_radius;
  result.control_count = std::max(4, std::min<int>(6, static_cast<int>(distance.size())));
  result.knots = open_uniform_knots(result.control_count, result.degree);
  Eigen::MatrixXd a(static_cast<Eigen::Index>(distance.size() + 2), result.control_count);
  Eigen::VectorXd y(static_cast<Eigen::Index>(distance.size() + 2));
  for (std::size_t row = 0; row < order.size(); ++row) {
    const double t = std::clamp(distance[order[row]] / result.max_distance, 0.0, 1.0);
    for (int col = 0; col < result.control_count; ++col)
      a(static_cast<Eigen::Index>(row), col) = bspline_basis(col, result.degree, t,
                                                            result.knots, result.control_count);
    y(static_cast<Eigen::Index>(row)) = radius[order[row]];
  }
  // Give the synthetic zero-radius tip enough weight to remain closed and smooth.
  for (int extra = 0; extra < 2; ++extra) {
    const Eigen::Index row = static_cast<Eigen::Index>(distance.size() + extra);
    for (int col = 0; col < result.control_count; ++col)
      a(row, col) = 4.0 * bspline_basis(col, result.degree, 0.0,
                                       result.knots, result.control_count);
    y(row) = 0;
  }
  result.controls = a.colPivHouseholderQr().solve(y);
  result.controls[0] = 0;
  return result;
}

double evaluate_profile(const SplineProfile& profile, const double distance) {
  const double t = std::clamp(distance / profile.max_distance, 0.0, 1.0);
  double radius = 0;
  for (int i = 0; i < profile.control_count; ++i)
    radius += profile.controls[i] *
              bspline_basis(i, profile.degree, t, profile.knots, profile.control_count);
  // Prevent sparse/noisy endpoint samples from creating an inverted tip.
  return std::clamp(radius, 0.0, profile.boundary_radius);
}

double angle_degrees(const Vec3& a, const Vec3& b) {
  if (a.squaredNorm() < 1e-18 || b.squaredNorm() < 1e-18) return 0;
  const double cosine = std::clamp(a.normalized().dot(b.normalized()), -1.0, 1.0);
  return std::acos(cosine) * 180.0 / kPi;
}

struct EndpointCenterlineFit {
  Vec3 boundary = Vec3::Zero();
  Vec3 outward_tangent = Vec3::UnitZ();
  Vec3 curvature = Vec3::Zero();
  double parameter_end = 0;
  double target_arc_length = 0;
  std::vector<std::size_t> source_ring_indices;
};

std::vector<std::size_t> collect_endpoint_centerline_fit_rings(
    const std::vector<RingReliabilityMetrics>& metrics, const std::size_t boundary_index,
    const bool bottom, const std::size_t requested_count) {
  std::vector<std::size_t> indices;
  indices.push_back(boundary_index);
  for (std::size_t step = 1; indices.size() < requested_count; ++step) {
    if (bottom) {
      if (boundary_index + step >= metrics.size()) break;
      const std::size_t index = boundary_index + step;
      if (metrics[index].reliable_for_mesh_boundary) indices.push_back(index);
    } else {
      if (step > boundary_index) break;
      const std::size_t index = boundary_index - step;
      if (metrics[index].reliable_for_mesh_boundary) indices.push_back(index);
    }
  }
  if (indices.size() < 2) throw std::runtime_error("Too few reliable Rings for endpoint centerline fitting");
  return indices;
}

Vec3 evaluate_endpoint_centerline_raw(const EndpointCenterlineFit& fit, const double parameter) {
  return fit.boundary + parameter * fit.outward_tangent +
         parameter * parameter * fit.curvature;
}

Vec3 evaluate_endpoint_centerline_tangent(const EndpointCenterlineFit& fit,
                                          const double parameter) {
  return (fit.outward_tangent + 2.0 * parameter * fit.curvature).normalized();
}

double endpoint_curve_length(const EndpointCenterlineFit& fit, const double parameter_end,
                             const std::size_t segments = 64) {
  Vec3 previous = fit.boundary;
  double length = 0;
  for (std::size_t i = 1; i <= segments; ++i) {
    const double parameter = parameter_end * static_cast<double>(i) /
                             static_cast<double>(segments);
    const Vec3 current = evaluate_endpoint_centerline_raw(fit, parameter);
    length += (current - previous).norm();
    previous = current;
  }
  return length;
}

EndpointCenterlineFit fit_endpoint_centerline(
    const std::vector<Slice>& rings, const std::vector<RingReliabilityMetrics>& metrics,
    const std::size_t boundary_index, const bool bottom,
    const RingReconstructionConfig& config, const double target_length) {
  EndpointCenterlineFit fit;
  fit.boundary = rings[boundary_index].center;
  fit.target_arc_length = target_length;
  fit.source_ring_indices = collect_endpoint_centerline_fit_rings(
      metrics, boundary_index, bottom, config.endpoint_centerline_fit_rings);
  std::vector<double> distances(fit.source_ring_indices.size(), 0.0);
  for (std::size_t i = 1; i < fit.source_ring_indices.size(); ++i)
    distances[i] = distances[i - 1] +
        (rings[fit.source_ring_indices[i]].center -
         rings[fit.source_ring_indices[i - 1]].center).norm();
  double s2 = 0, s3 = 0, s4 = 0;
  Vec3 rhs1 = Vec3::Zero(), rhs2 = Vec3::Zero();
  for (std::size_t i = 1; i < fit.source_ring_indices.size(); ++i) {
    const double s = distances[i];
    const Vec3 delta = rings[fit.source_ring_indices[i]].center - fit.boundary;
    s2 += s * s;
    s3 += s * s * s;
    s4 += s * s * s * s;
    rhs1 += s * delta;
    rhs2 += s * s * delta;
  }
  const double determinant = s2 * s4 - s3 * s3;
  Vec3 inward_tangent;
  Vec3 quadratic = Vec3::Zero();
  if (std::abs(determinant) > 1e-12) {
    inward_tangent = (s4 * rhs1 - s3 * rhs2) / determinant;
    quadratic = (-s3 * rhs1 + s2 * rhs2) / determinant;
  } else {
    inward_tangent = rings[fit.source_ring_indices.back()].center - fit.boundary;
  }
  if (inward_tangent.squaredNorm() < 1e-12)
    inward_tangent = bottom ? rings[boundary_index].tangent : -rings[boundary_index].tangent;
  inward_tangent.normalize();
  fit.outward_tangent = -inward_tangent;
  // s=-l during outward extrapolation, so the fitted quadratic term keeps its sign.
  fit.curvature = quadratic - fit.outward_tangent * quadratic.dot(fit.outward_tangent);
  const double allowed_turn = std::min(10.0, 0.5 * config.max_endpoint_turn_angle_deg);
  const double maximum_curvature = std::tan(allowed_turn * kPi / 180.0) /
                                   std::max(1e-9, 2.0 * target_length);
  if (fit.curvature.norm() > maximum_curvature)
    fit.curvature *= maximum_curvature / fit.curvature.norm();
  double low = 0, high = 1.5 * target_length;
  for (int iteration = 0; iteration < 32; ++iteration) {
    const double mid = 0.5 * (low + high);
    if (endpoint_curve_length(fit, mid) < target_length) low = mid;
    else high = mid;
  }
  fit.parameter_end = 0.5 * (low + high);
  return fit;
}

struct EndpointRadiusFit {
  double boundary_radius = 0;
  double outward_slope_times_length = 0;
};

EndpointRadiusFit fit_endpoint_radius_profile(
    const std::vector<Slice>& rings, const std::vector<std::size_t>& ring_indices,
    const double tip_length) {
  EndpointRadiusFit fit;
  fit.boundary_radius = rings[ring_indices.front()].representative_radius;
  double numerator = 0, denominator = 0, distance = 0;
  for (std::size_t i = 1; i < ring_indices.size(); ++i) {
    distance += (rings[ring_indices[i]].center - rings[ring_indices[i - 1]].center).norm();
    numerator += distance *
        (rings[ring_indices[i]].representative_radius - fit.boundary_radius);
    denominator += distance * distance;
  }
  const double inward_slope = denominator > 1e-12 ? numerator / denominator : 0;
  const double outward_slope = -std::max(0.0, inward_slope);
  fit.outward_slope_times_length = std::clamp(
      outward_slope * tip_length, -3.0 * fit.boundary_radius, 0.0);
  return fit;
}

double evaluate_endpoint_radius_profile(const EndpointRadiusFit& fit, const double t) {
  const double clamped = std::clamp(t, 0.0, 1.0);
  const double h00 = 2.0 * clamped * clamped * clamped -
                     3.0 * clamped * clamped + 1.0;
  const double h10 = clamped * clamped * clamped -
                     2.0 * clamped * clamped + clamped;
  return std::max(0.0, h00 * fit.boundary_radius +
                       h10 * fit.outward_slope_times_length);
}

using VertexRing = std::vector<Mesh::Vertex_index>;

VertexRing add_ring_vertices(Mesh& mesh, const Slice& ring) {
  VertexRing vertices;
  vertices.reserve(ring.radii.size());
  for (std::size_t j = 0; j < ring.radii.size(); ++j) {
    const double angle = 2.0 * kPi * static_cast<double>(j) /
                         static_cast<double>(ring.radii.size());
    const Vec3 p = ring.center + ring.radii[j] *
        (std::cos(angle) * ring.u + std::sin(angle) * ring.v);
    vertices.push_back(mesh.add_vertex(as_point(p)));
  }
  return vertices;
}

void connect_rings(Mesh& mesh, const VertexRing& lower, const VertexRing& upper) {
  const std::size_t n = lower.size();
  for (std::size_t j = 0; j < n; ++j) {
    const std::size_t next = (j + 1) % n;
    if (mesh.add_face(lower[j], lower[next], upper[next]) == Mesh::null_face() ||
        mesh.add_face(lower[j], upper[next], upper[j]) == Mesh::null_face())
      throw std::runtime_error("Failed to triangulate adjacent rings");
  }
}

}  // namespace

RingReconstructionResult reconstruct_ring_mesh(const std::vector<Point>& points,
                                               const RingReconstructionConfig& config) {
  if (points.size() < 100) throw std::runtime_error("Too few points for ring reconstruction");
  if (config.target_ring_count < 4 || config.angular_bins < 8)
    throw std::runtime_error("Ring count or angular bin count is too small");
  RingReconstructionResult output;
  output.target_ring_count = config.target_ring_count;
  auto& timing = output.timings;

  const auto pca_begin = Clock::now();
  const PcaFrame pca = compute_pca(points);
  const auto pca_end = Clock::now();
  timing.pca_ms = elapsed_ms(pca_begin, pca_end);

  const auto bin_begin = Clock::now();
  std::vector<double> projections(points.size());
  double s_min = std::numeric_limits<double>::infinity();
  double s_max = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < points.size(); ++i) {
    projections[i] = pca.axis.dot(as_vec(points[i]) - pca.centroid);
    s_min = std::min(s_min, projections[i]);
    s_max = std::max(s_max, projections[i]);
  }
  output.pca_length_mm = s_max - s_min;
  output.axial_step_mm = output.pca_length_mm / static_cast<double>(config.target_ring_count);
  std::vector<std::vector<std::size_t>> bins(config.target_ring_count);
  for (std::size_t i = 0; i < points.size(); ++i) {
    std::size_t bin = static_cast<std::size_t>(
        std::floor((projections[i] - s_min) / output.axial_step_mm));
    bin = std::min(bin, config.target_ring_count - 1);
    bins[bin].push_back(i);
  }
  std::vector<Slice> rings;
  rings.reserve(config.target_ring_count);
  const Vec3 pca_u = pca.transverse;
  const Vec3 pca_v = pca.axis.cross(pca_u).normalized();
  for (std::size_t b = 0; b < bins.size(); ++b) {
    if (bins[b].size() < config.min_slice_points) continue;
    Slice ring;
    ring.s = s_min + (static_cast<double>(b) + 0.5) * output.axial_step_mm;
    ring.point_indices = std::move(bins[b]);
    // Algebraic circle center in the PCA transverse plane is considerably less
    // sensitive to non-uniform scan density than the arithmetic point centroid.
    Eigen::Matrix3d normal = Eigen::Matrix3d::Zero();
    Vec3 rhs = Vec3::Zero();
    for (const auto index : ring.point_indices) {
      const Vec3 d = as_vec(points[index]) - pca.centroid;
      const double x = d.dot(pca_u);
      const double y = d.dot(pca_v);
      const Vec3 row(2.0 * x, 2.0 * y, 1.0);
      normal.noalias() += row * row.transpose();
      rhs.noalias() += row * (x * x + y * y);
    }
    const Vec3 circle = normal.ldlt().solve(rhs);
    if (circle.allFinite())
      ring.raw_center = pca.centroid + ring.s * pca.axis + circle.x() * pca_u + circle.y() * pca_v;
    else {
      for (const auto index : ring.point_indices) ring.raw_center += as_vec(points[index]);
      ring.raw_center /= static_cast<double>(ring.point_indices.size());
    }
    rings.push_back(std::move(ring));
  }
  const auto bin_end = Clock::now();
  timing.axial_binning_ms = elapsed_ms(bin_begin, bin_end);
  if (rings.size() < 4) throw std::runtime_error("Too few valid axial rings");
  output.actual_ring_count = rings.size();

  const auto center_begin = Clock::now();
  constexpr int weights[] = {1, 2, 3, 2, 1};
  for (std::size_t i = 0; i < rings.size(); ++i) {
    Vec3 sum = Vec3::Zero();
    int weight_sum = 0;
    for (int d = -2; d <= 2; ++d) {
      const long long j = static_cast<long long>(i) + d;
      if (j < 0 || j >= static_cast<long long>(rings.size())) continue;
      const int w = weights[d + 2];
      sum += static_cast<double>(w) * rings[static_cast<std::size_t>(j)].raw_center;
      weight_sum += w;
    }
    rings[i].center = sum / static_cast<double>(weight_sum);
  }
  const auto center_end = Clock::now();
  timing.centerline_ms = elapsed_ms(center_begin, center_end);

  const auto frame_begin = Clock::now();
  Vec3 previous_u = pca.transverse;
  for (std::size_t i = 0; i < rings.size(); ++i) {
    Vec3 tangent;
    if (i == 0) tangent = rings[1].center - rings[0].center;
    else if (i + 1 == rings.size()) tangent = rings[i].center - rings[i - 1].center;
    else tangent = rings[i + 1].center - rings[i - 1].center;
    if (tangent.dot(pca.axis) < 0) tangent = -tangent;
    tangent.normalize();
    Vec3 u = previous_u - tangent * previous_u.dot(tangent);
    if (u.squaredNorm() < 1e-12) u = tangent.unitOrthogonal();
    u.normalize();
    Vec3 v = tangent.cross(u).normalized();
    // u x v = tangent, so increasing sector angle is consistently oriented.
    rings[i].tangent = tangent;
    rings[i].u = u;
    rings[i].v = v;
    previous_u = u;
  }
  const auto frame_end = Clock::now();
  timing.local_frame_ms = elapsed_ms(frame_begin, frame_end);

  const auto radial_begin = Clock::now();
  std::vector<Slice> radial_valid;
  radial_valid.reserve(rings.size());
  for (auto& ring : rings) {
    std::vector<std::vector<double>> sectors(config.angular_bins);
    for (const auto index : ring.point_indices) {
      const Vec3 d = as_vec(points[index]) - ring.center;
      const double x = d.dot(ring.u);
      const double y = d.dot(ring.v);
      const double radius = std::hypot(x, y);
      double angle = std::atan2(y, x);
      if (angle < 0) angle += 2.0 * kPi;
      const std::size_t sector = std::min(config.angular_bins - 1,
          static_cast<std::size_t>(angle / (2.0 * kPi) * config.angular_bins));
      sectors[sector].push_back(radius);
    }
    ring.original_radii.assign(config.angular_bins, std::numeric_limits<double>::quiet_NaN());
    ring.radii.assign(config.angular_bins, std::numeric_limits<double>::quiet_NaN());
    std::size_t populated = 0;
    const auto mad_begin = Clock::now();
    for (std::size_t sector = 0; sector < sectors.size(); ++sector) {
      auto& values = sectors[sector];
      if (values.empty()) continue;
      std::sort(values.begin(), values.end());
      const double original_inner = quantile_sorted(values, config.inner_quantile);
      const double original_outer = quantile_sorted(values, config.outer_quantile);
      ring.original_radii[sector] = 0.5 * (original_inner + original_outer);
      const std::vector<double> filtered = filter_sector_radial_outliers(
          values, config.sector_mad_multiplier);
      const double inner = quantile_sorted(filtered, config.inner_quantile);
      const double outer = quantile_sorted(filtered, config.outer_quantile);
      ring.radii[sector] = 0.5 * (inner + outer);
      ++populated;
    }
    const auto mad_end = Clock::now();
    if (config.sector_mad_multiplier > 0)
      timing.sector_mad_filter_ms += elapsed_ms(mad_begin, mad_end);
    if (populated < std::max<std::size_t>(3, config.angular_bins / 2)) continue;
    ring.populated_sectors = populated;
    interpolate_empty_sectors(ring.original_radii);
    interpolate_empty_sectors(ring.radii);
    if (config.ring_smoothing_window > 1) {
      const auto smooth_begin = Clock::now();
      ring.radii = smooth_ring_radius_circular(ring.radii, config.ring_smoothing_window);
      const auto smooth_end = Clock::now();
      timing.ring_smoothing_ms += elapsed_ms(smooth_begin, smooth_end);
    }
    std::vector<double> representative = ring.radii;
    std::sort(representative.begin(), representative.end());
    ring.representative_radius = config.use_robust_endpoint_shape
        ? std::accumulate(ring.radii.begin(), ring.radii.end(), 0.0) /
          static_cast<double>(ring.radii.size())
        : quantile_sorted(representative, 0.5);
    radial_valid.push_back(std::move(ring));
  }
  rings = std::move(radial_valid);
  output.actual_ring_count = rings.size();
  const auto radial_end = Clock::now();
  timing.radial_sector_ms = elapsed_ms(radial_begin, radial_end);
  if (rings.size() < 4) throw std::runtime_error("Too few rings have adequate angular coverage");

  const auto quality_begin = Clock::now();
  output.ring_metrics.reserve(rings.size());
  for (std::size_t i = 0; i < rings.size(); ++i)
    output.ring_metrics.push_back(compute_ring_reliability(
        rings[i].radii, i, rings[i].point_indices.size(),
        static_cast<double>(rings[i].populated_sectors) /
        static_cast<double>(config.angular_bins)));
  const auto quality_end = Clock::now();
  timing.ring_quality_ms = elapsed_ms(quality_begin, quality_end);

  const auto debug_vertices = [&](const Slice& ring, const std::vector<double>& radii) {
    std::vector<Point> vertices;
    vertices.reserve(radii.size());
    for (std::size_t j = 0; j < radii.size(); ++j) {
      const double angle = 2.0 * kPi * static_cast<double>(j) /
                           static_cast<double>(radii.size());
      vertices.push_back(as_point(ring.center + radii[j] *
          (std::cos(angle) * ring.u + std::sin(angle) * ring.v)));
    }
    return vertices;
  };
  const std::size_t debug_count = std::min(config.endpoint_candidate_rings, rings.size());
  for (std::size_t i = 0; i < debug_count; ++i) {
    EndpointRingDebug bottom_debug;
    bottom_debug.ring_index = i;
    bottom_debug.original_vertices = debug_vertices(rings[i], rings[i].original_radii);
    bottom_debug.cleaned_vertices = debug_vertices(rings[i], rings[i].radii);
    output.bottom_endpoint_debug.push_back(std::move(bottom_debug));
    const std::size_t top_index = rings.size() - 1 - i;
    EndpointRingDebug top_debug;
    top_debug.ring_index = top_index;
    top_debug.original_vertices = debug_vertices(rings[top_index], rings[top_index].original_radii);
    top_debug.cleaned_vertices = debug_vertices(rings[top_index], rings[top_index].radii);
    output.top_endpoint_debug.push_back(std::move(top_debug));
  }

  EndpointSelection bottom_selection, top_selection;
  std::vector<std::size_t> bottom_radius_indices, top_radius_indices;
  std::size_t bottom_boundary_index = 0;
  std::size_t top_boundary_index = rings.size() - 1;
  if (config.endpoint == EndpointCompletion::cubic_bspline &&
      config.endpoint_generated_rings > 0) {
    if (config.use_robust_endpoint_shape) {
      const auto boundary_begin = Clock::now();
      bottom_selection = select_endpoint_reliable_rings(
          rings, output.ring_metrics, true, config);
      top_selection = select_endpoint_reliable_rings(
          rings, output.ring_metrics, false, config);
      const auto boundary_end = Clock::now();
      timing.endpoint_boundary_selection_ms = elapsed_ms(boundary_begin, boundary_end);
      timing.shape_template_ms = timing.endpoint_boundary_selection_ms;
      output.bottom_shape_template = bottom_selection.shape_template;
      output.top_shape_template = top_selection.shape_template;
      if (config.use_reliable_endpoint_boundary) {
        bottom_boundary_index = bottom_selection.boundary_index;
        top_boundary_index = top_selection.boundary_index;
        for (std::size_t index = bottom_boundary_index;
             index <= top_boundary_index &&
             bottom_radius_indices.size() < config.endpoint_fit_rings; ++index)
          if (output.ring_metrics[index].reliable_for_radius_fit)
            bottom_radius_indices.push_back(index);
        for (std::size_t index = top_boundary_index + 1;
             index-- > bottom_boundary_index &&
             top_radius_indices.size() < config.endpoint_fit_rings;)
          if (output.ring_metrics[index].reliable_for_radius_fit)
            top_radius_indices.push_back(index);
      } else {
        bottom_radius_indices = bottom_selection.radius_indices;
        top_radius_indices = top_selection.radius_indices;
      }
    } else {
      const std::size_t k = std::min(config.endpoint_fit_rings, rings.size());
      for (std::size_t i = 0; i < k; ++i) {
        bottom_radius_indices.push_back(i);
        top_radius_indices.push_back(rings.size() - 1 - i);
      }
    }
  }
  if (bottom_boundary_index >= top_boundary_index)
    throw std::runtime_error("Reliable endpoint boundaries leave no body Rings");

  std::vector<VertexRing> base_vertices;
  base_vertices.reserve(top_boundary_index - bottom_boundary_index + 1);
  const auto vertex_begin = Clock::now();
  for (std::size_t index = bottom_boundary_index; index <= top_boundary_index; ++index)
    base_vertices.push_back(add_ring_vertices(output.mesh, rings[index]));
  const auto vertex_end = Clock::now();
  timing.ring_vertex_generation_ms = elapsed_ms(vertex_begin, vertex_end);

  const auto mesh_begin = Clock::now();
  for (std::size_t i = 1; i < base_vertices.size(); ++i)
    connect_rings(output.mesh, base_vertices[i - 1], base_vertices[i]);
  const auto mesh_end = Clock::now();
  timing.ring_mesh_connection_ms = elapsed_ms(mesh_begin, mesh_end);

  const bool centerline_fixed = config.use_reliable_endpoint_boundary &&
      config.use_robust_endpoint_shape &&
      config.endpoint == EndpointCompletion::cubic_bspline &&
      config.endpoint_generated_rings > 0;
  SplineProfile bottom_profile, top_profile;
  EndpointCenterlineFit bottom_centerline_fit, top_centerline_fit;
  EndpointRadiusFit bottom_radius_fit, top_radius_fit;
  if (config.endpoint == EndpointCompletion::cubic_bspline &&
      config.endpoint_generated_rings > 0) {
    if (centerline_fixed) {
      const auto centerline_fit_begin = Clock::now();
      const double bottom_length = config.tip_length_ratio *
          rings[bottom_boundary_index].representative_radius;
      const double top_length = config.tip_length_ratio *
          rings[top_boundary_index].representative_radius;
      bottom_centerline_fit = fit_endpoint_centerline(rings, output.ring_metrics,
          bottom_boundary_index, true, config, bottom_length);
      top_centerline_fit = fit_endpoint_centerline(rings, output.ring_metrics,
          top_boundary_index, false, config, top_length);
      const auto centerline_fit_end = Clock::now();
      timing.endpoint_centerline_fit_ms = elapsed_ms(
          centerline_fit_begin, centerline_fit_end);
      const auto radius_fit_begin = Clock::now();
      bottom_radius_fit = fit_endpoint_radius_profile(
          rings, bottom_radius_indices, bottom_length);
      top_radius_fit = fit_endpoint_radius_profile(
          rings, top_radius_indices, top_length);
      const auto radius_fit_end = Clock::now();
      timing.endpoint_radius_fit_ms = elapsed_ms(radius_fit_begin, radius_fit_end);
      timing.bottom_spline_fit_ms = 0.5 * timing.endpoint_radius_fit_ms;
      timing.top_spline_fit_ms = 0.5 * timing.endpoint_radius_fit_ms;
    } else {
      const auto bottom_fit_begin = Clock::now();
      bottom_profile = fit_endpoint(rings, true, bottom_radius_indices, s_min,
          config.use_robust_endpoint_shape && config.enforce_endpoint_monotonic_radius,
          config.use_robust_endpoint_shape ? &output.ring_metrics : nullptr);
      const auto bottom_fit_end = Clock::now();
      timing.bottom_spline_fit_ms = elapsed_ms(bottom_fit_begin, bottom_fit_end);
      const auto top_fit_begin = Clock::now();
      top_profile = fit_endpoint(rings, false, top_radius_indices, s_max,
          config.use_robust_endpoint_shape && config.enforce_endpoint_monotonic_radius,
          config.use_robust_endpoint_shape ? &output.ring_metrics : nullptr);
      const auto top_fit_end = Clock::now();
      timing.top_spline_fit_ms = elapsed_ms(top_fit_begin, top_fit_end);
      timing.endpoint_radius_fit_ms = timing.bottom_spline_fit_ms + timing.top_spline_fit_ms;
    }
  }

  std::vector<double> bottom_profile_values, top_profile_values;
  std::vector<Vec3> bottom_centers, top_centers;
  std::vector<Vec3> bottom_outward_tangents, top_outward_tangents;
  Vec3 fitted_bottom_tip = Vec3::Zero(), fitted_top_tip = Vec3::Zero();
  if (config.endpoint == EndpointCompletion::cubic_bspline &&
      config.endpoint_generated_rings > 0) {
    bottom_profile_values.resize(config.endpoint_generated_rings);
    top_profile_values.resize(config.endpoint_generated_rings);
    if (centerline_fixed) {
      const auto centerline_eval_begin = Clock::now();
      bottom_centers.resize(config.endpoint_generated_rings);
      top_centers.resize(config.endpoint_generated_rings);
      bottom_outward_tangents.resize(config.endpoint_generated_rings);
      top_outward_tangents.resize(config.endpoint_generated_rings);
      double previous_bottom_radius = bottom_radius_fit.boundary_radius;
      double previous_top_radius = top_radius_fit.boundary_radius;
      for (std::size_t i = 0; i < config.endpoint_generated_rings; ++i) {
        const double t = static_cast<double>(i + 1) /
                         static_cast<double>(config.endpoint_generated_rings + 1);
        const double bottom_parameter = t * bottom_centerline_fit.parameter_end;
        const double top_parameter = t * top_centerline_fit.parameter_end;
        bottom_centers[i] = evaluate_endpoint_centerline_raw(
            bottom_centerline_fit, bottom_parameter);
        top_centers[i] = evaluate_endpoint_centerline_raw(
            top_centerline_fit, top_parameter);
        bottom_outward_tangents[i] = evaluate_endpoint_centerline_tangent(
            bottom_centerline_fit, bottom_parameter);
        top_outward_tangents[i] = evaluate_endpoint_centerline_tangent(
            top_centerline_fit, top_parameter);
        double bottom_radius = evaluate_endpoint_radius_profile(bottom_radius_fit, t);
        double top_radius = evaluate_endpoint_radius_profile(top_radius_fit, t);
        if (config.enforce_endpoint_monotonic_radius) {
          bottom_radius = std::min(bottom_radius, previous_bottom_radius);
          top_radius = std::min(top_radius, previous_top_radius);
        }
        bottom_profile_values[i] = std::max(0.0, bottom_radius);
        top_profile_values[i] = std::max(0.0, top_radius);
        previous_bottom_radius = bottom_profile_values[i];
        previous_top_radius = top_profile_values[i];
      }
      fitted_bottom_tip = evaluate_endpoint_centerline_raw(
          bottom_centerline_fit, bottom_centerline_fit.parameter_end);
      fitted_top_tip = evaluate_endpoint_centerline_raw(
          top_centerline_fit, top_centerline_fit.parameter_end);
      const auto centerline_eval_end = Clock::now();
      timing.endpoint_centerline_eval_ms = elapsed_ms(
          centerline_eval_begin, centerline_eval_end);
    } else {
      const auto bottom_eval_begin = Clock::now();
      for (std::size_t i = 0; i < bottom_profile_values.size(); ++i) {
        const double fraction = static_cast<double>(i + 1) /
                                static_cast<double>(bottom_profile_values.size() + 1);
        const double radius = evaluate_profile(bottom_profile,
            fraction * bottom_profile.boundary_distance);
        bottom_profile_values[i] = config.use_robust_endpoint_shape ? radius
            : radius / std::max(1e-9, bottom_profile.boundary_radius);
      }
      const auto bottom_eval_end = Clock::now();
      timing.bottom_spline_eval_ms = elapsed_ms(bottom_eval_begin, bottom_eval_end);
      const auto top_eval_begin = Clock::now();
      for (std::size_t i = 0; i < top_profile_values.size(); ++i) {
        const double fraction = static_cast<double>(i + 1) /
                                static_cast<double>(top_profile_values.size() + 1);
        const double radius = evaluate_profile(top_profile,
            fraction * top_profile.boundary_distance);
        top_profile_values[i] = config.use_robust_endpoint_shape ? radius
            : radius / std::max(1e-9, top_profile.boundary_radius);
      }
      const auto top_eval_end = Clock::now();
      timing.top_spline_eval_ms = elapsed_ms(top_eval_begin, top_eval_end);
    }
  }

  const auto endpoint_gen_begin = Clock::now();
  const Slice& bottom = rings[bottom_boundary_index];
  const Slice& top = rings[top_boundary_index];
  const double bottom_gap = std::max(0.0, bottom.s - s_min);
  const double top_gap = std::max(0.0, s_max - top.s);
  const bool spline_endpoint = config.endpoint == EndpointCompletion::cubic_bspline &&
                               config.endpoint_generated_rings > 0;
  // The baseline is a true coplanar triangle fan. Only the spline variant
  // extends generated rings from the last reliable ring toward the PCA extreme.
  const Vec3 bottom_tip_position = centerline_fixed ? fitted_bottom_tip :
      (spline_endpoint ? bottom.center - bottom.tangent * bottom_gap : bottom.center);
  const Vec3 top_tip_position = centerline_fixed ? fitted_top_tip :
      (spline_endpoint ? top.center + top.tangent * top_gap : top.center);
  std::vector<VertexRing> bottom_generated, top_generated;
  if (!bottom_profile_values.empty()) {
    if (centerline_fixed) {
      const auto set_frame = [](Slice& generated, const Slice& boundary,
                                const Vec3& increasing_tangent) {
        generated.tangent = increasing_tangent.normalized();
        generated.u = boundary.u - generated.tangent * boundary.u.dot(generated.tangent);
        if (generated.u.squaredNorm() < 1e-12)
          generated.u = generated.tangent.unitOrthogonal();
        generated.u.normalize();
        generated.v = generated.tangent.cross(generated.u).normalized();
      };
      std::vector<Slice> bottom_slices;
      bottom_slices.reserve(config.endpoint_generated_rings);
      for (std::size_t i = 0; i < config.endpoint_generated_rings; ++i) {
        const double t = static_cast<double>(i + 1) /
                         static_cast<double>(config.endpoint_generated_rings + 1);
        Slice generated = bottom;
        generated.center = bottom_centers[i];
        set_frame(generated, bottom, -bottom_outward_tangents[i]);
        generated.radii.resize(config.angular_bins);
        const double regularization = config.tip_shape_regularization ? t : 0.0;
        for (std::size_t sector = 0; sector < config.angular_bins; ++sector) {
          const double shape = (1.0 - regularization) *
              output.bottom_shape_template.normalized_radii[sector] + regularization;
          generated.radii[sector] = bottom_profile_values[i] * shape;
        }
        output.bottom_generated_metrics.push_back(compute_ring_reliability(
            generated.radii, i, 0, 1.0));
        bottom_slices.push_back(std::move(generated));
      }
      for (auto iterator = bottom_slices.rbegin(); iterator != bottom_slices.rend(); ++iterator)
        bottom_generated.push_back(add_ring_vertices(output.mesh, *iterator));
      for (std::size_t i = 0; i < config.endpoint_generated_rings; ++i) {
        const double t = static_cast<double>(i + 1) /
                         static_cast<double>(config.endpoint_generated_rings + 1);
        Slice generated = top;
        generated.center = top_centers[i];
        set_frame(generated, top, top_outward_tangents[i]);
        generated.radii.resize(config.angular_bins);
        const double regularization = config.tip_shape_regularization ? t : 0.0;
        for (std::size_t sector = 0; sector < config.angular_bins; ++sector) {
          const double shape = (1.0 - regularization) *
              output.top_shape_template.normalized_radii[sector] + regularization;
          generated.radii[sector] = top_profile_values[i] * shape;
        }
        output.top_generated_metrics.push_back(compute_ring_reliability(
            generated.radii, i, 0, 1.0));
        top_generated.push_back(add_ring_vertices(output.mesh, generated));
      }
    } else {
      for (std::size_t i = 0; i < bottom_profile_values.size(); ++i) {
        const double fraction = static_cast<double>(i + 1) /
                                static_cast<double>(bottom_profile_values.size() + 1);
        Slice generated = bottom;
        generated.center = bottom_tip_position + fraction * (bottom.center - bottom_tip_position);
        if (config.use_robust_endpoint_shape) {
          generated.radii.resize(config.angular_bins);
          const double regularization = config.tip_shape_regularization ? 1.0 - fraction : 0.0;
          for (std::size_t sector = 0; sector < config.angular_bins; ++sector) {
            const double shape = (1.0 - regularization) *
                output.bottom_shape_template.normalized_radii[sector] + regularization;
            generated.radii[sector] = bottom_profile_values[i] * shape;
          }
        } else {
          for (double& radius : generated.radii) radius *= bottom_profile_values[i];
        }
        output.bottom_generated_metrics.push_back(compute_ring_reliability(
            generated.radii, i, 0, 1.0));
        bottom_generated.push_back(add_ring_vertices(output.mesh, generated));
      }
      for (std::size_t i = 0; i < top_profile_values.size(); ++i) {
        const double fraction = static_cast<double>(i + 1) /
                                static_cast<double>(top_profile_values.size() + 1);
        Slice generated = top;
        generated.center = top.center + fraction * (top_tip_position - top.center);
        const double profile_value = top_profile_values[top_profile_values.size() - 1 - i];
        if (config.use_robust_endpoint_shape) {
          generated.radii.resize(config.angular_bins);
          const double regularization = config.tip_shape_regularization ? fraction : 0.0;
          for (std::size_t sector = 0; sector < config.angular_bins; ++sector) {
            const double shape = (1.0 - regularization) *
                output.top_shape_template.normalized_radii[sector] + regularization;
            generated.radii[sector] = profile_value * shape;
          }
        } else {
          for (double& radius : generated.radii) radius *= profile_value;
        }
        output.top_generated_metrics.push_back(compute_ring_reliability(
            generated.radii, i, 0, 1.0));
        top_generated.push_back(add_ring_vertices(output.mesh, generated));
      }
    }
  }
  const auto bottom_tip = output.mesh.add_vertex(as_point(bottom_tip_position));
  const auto top_tip = output.mesh.add_vertex(as_point(top_tip_position));
  const auto endpoint_gen_end = Clock::now();
  timing.endpoint_ring_generation_ms = elapsed_ms(endpoint_gen_begin, endpoint_gen_end);
  timing.endpoint_generation_ms = timing.endpoint_ring_generation_ms;

  const auto endpoint_connect_begin = Clock::now();
  const VertexRing& first_bottom_ring = bottom_generated.empty() ? base_vertices.front()
                                                                 : bottom_generated.front();
  for (std::size_t j = 0; j < config.angular_bins; ++j) {
    const std::size_t next = (j + 1) % config.angular_bins;
    if (output.mesh.add_face(bottom_tip, first_bottom_ring[next], first_bottom_ring[j]) ==
        Mesh::null_face()) throw std::runtime_error("Failed to create bottom tip fan");
  }
  for (std::size_t i = 1; i < bottom_generated.size(); ++i)
    connect_rings(output.mesh, bottom_generated[i - 1], bottom_generated[i]);
  if (!bottom_generated.empty()) connect_rings(output.mesh, bottom_generated.back(), base_vertices.front());

  if (!top_generated.empty()) {
    connect_rings(output.mesh, base_vertices.back(), top_generated.front());
    for (std::size_t i = 1; i < top_generated.size(); ++i)
      connect_rings(output.mesh, top_generated[i - 1], top_generated[i]);
  }
  const VertexRing& last_top_ring = top_generated.empty() ? base_vertices.back()
                                                           : top_generated.back();
  for (std::size_t j = 0; j < config.angular_bins; ++j) {
    const std::size_t next = (j + 1) % config.angular_bins;
    if (output.mesh.add_face(last_top_ring[j], last_top_ring[next], top_tip) == Mesh::null_face())
      throw std::runtime_error("Failed to create top tip fan");
  }
  const auto endpoint_connect_end = Clock::now();
  timing.endpoint_mesh_connection_ms = elapsed_ms(endpoint_connect_begin, endpoint_connect_end);

  if (centerline_fixed) {
    const auto fill_geometry_metrics = [&](const bool is_bottom,
                                           const EndpointCenterlineFit& fit,
                                           const EndpointRadiusFit& radius_fit,
                                           const std::vector<Vec3>& centers,
                                           const Vec3& tip,
                                           const std::vector<double>& generated_radii,
                                           EndpointGeometryMetrics& metrics,
                                           std::vector<EndpointCenterlinePoint>& debug) {
      const std::size_t boundary_index = is_bottom ? bottom_boundary_index : top_boundary_index;
      metrics.boundary_ring_index = boundary_index;
      metrics.boundary_radius = radius_fit.boundary_radius;
      metrics.tip_length_ratio = config.tip_length_ratio;
      metrics.generated_ring_count = centers.size();
      const Vec3 first_segment = centers.front() - fit.boundary;
      const Vec3 first_tangent = first_segment.normalized();
      metrics.body_tangent = {{fit.outward_tangent.x(), fit.outward_tangent.y(),
                               fit.outward_tangent.z()}};
      metrics.first_generated_tangent = {{first_tangent.x(), first_tangent.y(),
                                          first_tangent.z()}};
      metrics.boundary_tangent_angle_deg = angle_degrees(fit.outward_tangent, first_segment);
      std::vector<Vec3> path;
      path.push_back(fit.boundary);
      path.insert(path.end(), centers.begin(), centers.end());
      path.push_back(tip);
      metrics.tip_length_mm = 0;
      for (std::size_t i = 1; i < path.size(); ++i)
        metrics.tip_length_mm += (path[i] - path[i - 1]).norm();
      metrics.max_endpoint_centerline_turn_deg = 0;
      for (std::size_t i = 1; i + 1 < path.size(); ++i)
        metrics.max_endpoint_centerline_turn_deg = std::max(
            metrics.max_endpoint_centerline_turn_deg,
            angle_degrees(path[i] - path[i - 1], path[i + 1] - path[i]));
      metrics.generated_radius_monotonic = true;
      double previous_radius = radius_fit.boundary_radius;
      for (const double radius : generated_radii) {
        if (radius > previous_radius + 1e-9) metrics.generated_radius_monotonic = false;
        previous_radius = radius;
      }
      if (previous_radius < -1e-9) metrics.generated_radius_monotonic = false;
      metrics.tangent_warning = metrics.boundary_tangent_angle_deg >
                                config.max_boundary_tangent_angle_deg;
      metrics.turn_warning = metrics.max_endpoint_centerline_turn_deg >
                             config.max_endpoint_turn_angle_deg;
      std::size_t debug_index = 0;
      for (auto iterator = fit.source_ring_indices.rbegin();
           iterator != fit.source_ring_indices.rend(); ++iterator) {
        if (*iterator == boundary_index) continue;
        debug.push_back({debug_index++, "body_ring", as_point(rings[*iterator].center),
                         rings[*iterator].representative_radius});
      }
      debug.push_back({debug_index++, "boundary", as_point(fit.boundary),
                       radius_fit.boundary_radius});
      for (std::size_t i = 0; i < centers.size(); ++i)
        debug.push_back({debug_index++, "generated_ring", as_point(centers[i]),
                         generated_radii[i]});
      debug.push_back({debug_index, "tip", as_point(tip), 0.0});
    };
    fill_geometry_metrics(true, bottom_centerline_fit, bottom_radius_fit,
        bottom_centers, bottom_tip_position, bottom_profile_values,
        output.bottom_geometry_metrics, output.bottom_endpoint_centerline);
    fill_geometry_metrics(false, top_centerline_fit, top_radius_fit,
        top_centers, top_tip_position, top_profile_values,
        output.top_geometry_metrics, output.top_endpoint_centerline);
  }

  timing.endpoint_completion_ms = timing.shape_template_ms +
      timing.endpoint_centerline_fit_ms + timing.endpoint_centerline_eval_ms +
      timing.bottom_spline_fit_ms + timing.top_spline_fit_ms +
      timing.bottom_spline_eval_ms + timing.top_spline_eval_ms +
      timing.endpoint_ring_generation_ms + timing.endpoint_mesh_connection_ms;
  timing.nurbs_endpoint_total_ms = timing.shape_template_ms +
      timing.endpoint_centerline_fit_ms + timing.endpoint_centerline_eval_ms +
      timing.bottom_spline_fit_ms + timing.top_spline_fit_ms +
      timing.bottom_spline_eval_ms + timing.top_spline_eval_ms +
      (config.endpoint == EndpointCompletion::cubic_bspline
           ? timing.endpoint_ring_generation_ms + timing.endpoint_mesh_connection_ms : 0.0);
  timing.ring_reconstruction_ms = timing.pca_ms + timing.axial_binning_ms +
      timing.centerline_ms + timing.local_frame_ms + timing.radial_sector_ms +
      timing.ring_quality_ms +
      timing.ring_vertex_generation_ms + timing.ring_mesh_connection_ms +
      timing.endpoint_completion_ms;
  return output;
}
