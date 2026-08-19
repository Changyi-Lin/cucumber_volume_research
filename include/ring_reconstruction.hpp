#pragma once

#include "types.hpp"

#include <cstddef>
#include <vector>

enum class EndpointCompletion { flat_cap, cubic_bspline };

struct RingReconstructionConfig {
  std::size_t target_ring_count = 90;
  std::size_t angular_bins = 48;
  std::size_t min_slice_points = 48;
  double inner_quantile = 0.10;
  double outer_quantile = 0.90;
  EndpointCompletion endpoint = EndpointCompletion::flat_cap;
  std::size_t endpoint_fit_rings = 8;
  std::size_t endpoint_generated_rings = 4;
  double sector_mad_multiplier = 3.0;
  std::size_t ring_smoothing_window = 5;
  std::size_t endpoint_candidate_rings = 10;
  std::size_t endpoint_shape_ring_count = 4;
  std::size_t endpoint_centerline_fit_rings = 5;
  double endpoint_absolute_max_cv = 0.12;
  double endpoint_absolute_max_normalized_roughness = 0.12;
  double endpoint_quality_mad_multiplier = 2.5;
  bool enforce_endpoint_monotonic_radius = true;
  bool tip_shape_regularization = false;
  bool use_robust_endpoint_shape = true;
  bool use_reliable_endpoint_boundary = true;
  double tip_length_ratio = 0.8;
  double max_boundary_tangent_angle_deg = 15.0;
  double max_endpoint_turn_angle_deg = 20.0;
};

struct RingReliabilityMetrics {
  std::size_t ring_index = 0;
  std::size_t point_count = 0;
  double angular_coverage = 0;
  double mean_radius = 0;
  double radius_stddev = 0;
  double radius_cv = 0;
  double angular_roughness = 0;
  double normalized_roughness = 0;
  bool reliable_for_shape = false;
  bool reliable_for_radius_fit = false;
  bool reliable_for_mesh_boundary = false;
  bool endpoint_radius_rebound = false;
};

struct EndpointShapeTemplate {
  std::vector<std::size_t> source_ring_indices;
  std::vector<double> normalized_radii;
};

struct EndpointRingDebug {
  std::size_t ring_index = 0;
  std::vector<Point> original_vertices;
  std::vector<Point> cleaned_vertices;
};

struct EndpointGeometryMetrics {
  std::size_t boundary_ring_index = 0;
  double boundary_radius = 0;
  double tip_length_mm = 0;
  double tip_length_ratio = 0;
  std::array<double, 3> body_tangent{{0, 0, 0}};
  std::array<double, 3> first_generated_tangent{{0, 0, 0}};
  double boundary_tangent_angle_deg = 0;
  double max_endpoint_centerline_turn_deg = 0;
  std::size_t generated_ring_count = 0;
  bool generated_radius_monotonic = false;
  bool tangent_warning = false;
  bool turn_warning = false;
};

struct EndpointCenterlinePoint {
  std::size_t index = 0;
  std::string type;
  Point center;
  double radius = 0;
};

struct RingReconstructionTimings {
  double pca_ms = 0;
  double axial_binning_ms = 0;
  double centerline_ms = 0;
  double local_frame_ms = 0;
  double radial_sector_ms = 0;
  double sector_mad_filter_ms = 0;
  double ring_smoothing_ms = 0;
  double ring_quality_ms = 0;
  double ring_vertex_generation_ms = 0;
  double ring_mesh_connection_ms = 0;
  double bottom_spline_fit_ms = 0;
  double top_spline_fit_ms = 0;
  double bottom_spline_eval_ms = 0;
  double top_spline_eval_ms = 0;
  double endpoint_ring_generation_ms = 0;
  double endpoint_generation_ms = 0;
  double endpoint_mesh_connection_ms = 0;
  double endpoint_completion_ms = 0;
  double shape_template_ms = 0;
  double endpoint_boundary_selection_ms = 0;
  double endpoint_centerline_fit_ms = 0;
  double endpoint_centerline_eval_ms = 0;
  double endpoint_radius_fit_ms = 0;
  double nurbs_endpoint_total_ms = 0;
  double ring_reconstruction_ms = 0;
};

struct RingReconstructionResult {
  Mesh mesh;
  RingReconstructionTimings timings;
  std::size_t target_ring_count = 0;
  std::size_t actual_ring_count = 0;
  double axial_step_mm = 0;
  double pca_length_mm = 0;
  std::vector<RingReliabilityMetrics> ring_metrics;
  EndpointShapeTemplate bottom_shape_template;
  EndpointShapeTemplate top_shape_template;
  std::vector<EndpointRingDebug> bottom_endpoint_debug;
  std::vector<EndpointRingDebug> top_endpoint_debug;
  std::vector<RingReliabilityMetrics> bottom_generated_metrics;
  std::vector<RingReliabilityMetrics> top_generated_metrics;
  EndpointGeometryMetrics bottom_geometry_metrics;
  EndpointGeometryMetrics top_geometry_metrics;
  std::vector<EndpointCenterlinePoint> bottom_endpoint_centerline;
  std::vector<EndpointCenterlinePoint> top_endpoint_centerline;
};

RingReconstructionResult reconstruct_ring_mesh(
    const std::vector<Point>& points, const RingReconstructionConfig& config);
