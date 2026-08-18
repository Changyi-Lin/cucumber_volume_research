#include "alpha_wrap.hpp"
#include "denoise.hpp"
#include "mesh_validation.hpp"
#include "ply_reader.hpp"
#include "volume.hpp"
#include "voxel_downsample.hpp"

#include <CGAL/IO/polygon_mesh_io.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

namespace {

double elapsed_ms(const Clock::time_point begin, const Clock::time_point end) {
  return std::chrono::duration<double, std::milli>(end - begin).count();
}

struct Arguments {
  std::unordered_map<std::string, std::string> values;
  std::vector<std::string> flags;

  bool has(const std::string& key) const {
    return values.count(key) || std::find(flags.begin(), flags.end(), key) != flags.end();
  }
  std::string get(const std::string& key, const std::string& fallback = {}) const {
    const auto it = values.find(key);
    return it == values.end() ? fallback : it->second;
  }
};

Arguments parse_arguments(const int argc, char** argv) {
  Arguments args;
  for (int i = 1; i < argc; ++i) {
    std::string key = argv[i];
    if (key.rfind("--", 0) != 0) throw std::runtime_error("Unexpected argument: " + key);
    if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0)
      args.values[key] = argv[++i];
    else
      args.flags.push_back(key);
  }
  return args;
}

std::vector<double> parse_doubles(const std::string& text) {
  std::vector<double> values;
  std::istringstream stream(text);
  std::string token;
  while (std::getline(stream, token, ',')) values.push_back(std::stod(token));
  return values;
}

fs::path find_default_input() {
  for (const auto& candidate : {fs::path("pointcloud_3100.ply"), fs::path("../pointcloud_3100.ply")})
    if (fs::exists(candidate)) return candidate;
  for (const auto& root : {fs::current_path(), fs::current_path().parent_path()}) {
    for (const auto& entry : fs::recursive_directory_iterator(root))
      if (entry.is_regular_file() && entry.path().filename() == "pointcloud_3100.ply") return entry.path();
  }
  throw std::runtime_error("pointcloud_3100.ply was not found");
}

std::string clean_csv(std::string value) {
  std::replace(value.begin(), value.end(), ',', ';');
  std::replace(value.begin(), value.end(), '\n', ' ');
  return value;
}

void print_stats(const DatasetStats& s, const double load_ms) {
  std::cout << std::fixed << std::setprecision(3)
            << "Dataset points: " << s.point_count << "\n"
            << "AABB min: " << s.min[0] << ", " << s.min[1] << ", " << s.min[2] << " mm\n"
            << "AABB max: " << s.max[0] << ", " << s.max[1] << ", " << s.max[2] << " mm\n"
            << "AABB diagonal: " << s.aabb_diagonal << " mm\n"
            << "PCA extents (ascending axes): " << s.pca_extents[0] << " x "
            << s.pca_extents[1] << " x " << s.pca_extents[2] << " mm\n"
            << "PLY load: " << load_ms << " ms\n";
}

struct RunOutput {
  BenchmarkResult result;
  Mesh mesh;
  DenoiseConfig denoise_config;
  DenoiseMetrics denoise_metrics;
  std::vector<Point> denoised_points;
  std::vector<Point> removed_points;
  std::vector<Point> voxel_points;
  DenoiseStages stages;
  std::vector<ComponentInfo> components;
  std::vector<std::size_t> component_labels;
  std::size_t points_before_denoise = 0;
  double injected_outlier_percent = 0;
  double injected_outlier_distance_mm = 0;
  std::string injected_outlier_mode = "dispersed";
};

std::string denoise_pipeline_name(const DenoiseConfig& config) {
  std::string name = "VOXEL";
  if (config.method == DenoiseMethod::sor) name += "+SOR";
  else if (config.method == DenoiseMethod::ror) name += "+ROR";
  else if (config.method == DenoiseMethod::sor_ror) name += "+SOR+ROR";
  if (config.clustering) name += "+COMPONENT";
  if (config.adaptive_local) name += "+ADAPTIVE_LOCAL";
  return name;
}

void write_result_csv(const fs::path& path, const BenchmarkResult& r) {
  fs::create_directories(path.parent_path().empty() ? fs::path(".") : path.parent_path());
  const bool needs_header = !fs::exists(path) || fs::file_size(path) == 0;
  std::ofstream out(path, std::ios::app);
  if (!out) throw std::runtime_error("Cannot write CSV: " + path.string());
  if (needs_header) {
    out << "voxel_size_mm,original_points,input_points,alpha_mm,offset_mm,run_index,"
           "ply_load_ms,downsample_ms,alpha_wrap_ms,component_cleanup_ms,validation_ms,volume_ms,total_ms,"
           "reconstruction_only_ms,online_core_ms,vertices,faces,edges,connected_components,"
           "watertight,is_triangle_mesh,self_intersection,manifold,outward_oriented,"
           "two_sided_wrap,euler_characteristic,genus,median_line_intersections,"
           "lines_with_four_or_more,tested_lines,signed_volume_mm3,volume_mm3,volume_ml,"
           "surface_area_mm2,raw_connected_components,removed_components,"
           "largest_component_face_ratio,valid,failure_reason\n";
  }
  const auto& v = r.validation;
  out << std::setprecision(12)
      << r.voxel_size_mm << ',' << r.original_points << ',' << r.input_points << ','
      << r.alpha_mm << ',' << r.offset_mm << ',' << r.run_index << ',' << r.ply_load_ms << ','
      << r.downsample_ms << ',' << r.alpha_wrap_ms << ',' << r.component_cleanup_ms << ',' << r.validation_ms << ','
      << r.volume_ms << ',' << r.total_ms << ',' << r.reconstruction_only_ms << ','
      << r.online_core_ms << ',' << v.vertices << ',' << v.faces << ',' << v.edges << ','
      << v.connected_components << ',' << v.is_closed << ',' << v.is_triangle_mesh << ','
      << v.self_intersection << ',' << v.manifold << ',' << v.outward_oriented << ','
      << v.two_sided_wrap << ',' << v.euler_characteristic << ',' << v.genus << ','
      << v.median_line_intersections << ',' << v.lines_with_four_or_more << ','
      << v.tested_lines << ',' << r.signed_volume_mm3 << ',' << r.volume_mm3 << ','
      << r.volume_mm3 / 1000.0 << ',' << v.surface_area_mm2 << ','
      << r.cleanup.raw_connected_components << ',' << r.cleanup.removed_components << ','
      << r.cleanup.largest_component_face_ratio << ',' << r.valid << ','
      << clean_csv(r.failure_reason) << '\n';
}

void write_denoise_result_csv(const fs::path& path, const RunOutput& output) {
  const BenchmarkResult& r = output.result;
  const DenoiseMetrics& d = output.denoise_metrics;
  const DenoiseConfig& c = output.denoise_config;
  fs::create_directories(path.parent_path().empty() ? fs::path(".") : path.parent_path());
  const bool needs_header = !fs::exists(path) || fs::file_size(path) == 0;
  std::ofstream out(path, std::ios::app);
  if (!out) throw std::runtime_error("Cannot write denoising CSV: " + path.string());
  if (needs_header) {
    out << "voxel_size_mm,voxel_size,denoise_method,denoise_pipeline,"
           "sor_k,sor_std_ratio,sor_std,ror_radius_mm,ror_radius,ror_min_neighbors,"
           "end_protection,end_region_fraction,clustering,cluster_radius_mm,component_radius,"
           "component_keep_largest,component_min_size,confirmed_disconnected_noise,"
           "adaptive_local,local_knn,local_pca_knn,"
           "local_mad_factor,local_regional_threshold,local_end_threshold_multiplier,local_score_threshold,"
           "injected_outlier_percent,injected_outlier_distance_mm,injected_outlier_mode,"
           "points_before,points_after,points_after_sor,points_after_ror,points_after_component,points_after_local,"
           "removed_points,removed_percent,removed_sor_points,removed_ror_points,"
           "removed_component_points,removed_local_points,"
           "top_before,top_after,top_removed_percent,middle_before,middle_after,middle_removed_percent,"
           "bottom_before,bottom_after,bottom_removed_percent,"
           "length_before_mm,length_after_mm,length_change_mm,length_change_percent,"
           "confirmed_component_length_change_mm,unexplained_length_change_mm,possible_over_denoise,"
           "warning_end_damage,invalid_over_filtering,"
           "number_of_clusters,component_count,largest_cluster_points,largest_component_points,"
           "second_cluster_points,removed_cluster_points,"
           "local_bottom_density_median,local_bottom_density_mad,local_bottom_density_threshold,"
           "local_middle_density_median,local_middle_density_mad,local_middle_density_threshold,"
           "local_top_density_median,local_top_density_mad,local_top_density_threshold,"
           "load_ms,voxel_ms,sor_ms,ror_ms,clustering_ms,component_ms,local_ms,denoise_total_ms,denoise_ms,"
           "alpha_mm,offset_mm,"
           "alpha_wrap_ms,component_cleanup_ms,validation_ms,volume_ms,deployment_total_ms,total_ms,research_total_ms,"
           "vertices,faces,connected_components,watertight,closed,manifold,self_intersection,"
           "outward_oriented,two_sided_wrap,"
           "volume_mm3,volume_ml,surface_area_mm2,valid,failure_reason\n";
  }
  const auto& v = r.validation;
  out << std::setprecision(12)
      << r.voxel_size_mm << ',' << r.voxel_size_mm << ',' << denoise_method_name(c.method) << ','
      << denoise_pipeline_name(c) << ',' << c.sor_k << ','
      << c.sor_std_ratio << ',' << c.sor_std_ratio << ',' << c.ror_radius_mm << ','
      << c.ror_radius_mm << ',' << c.ror_min_neighbors << ','
      << c.end_protection_fraction << ',' << c.metric_end_fraction << ',' << c.clustering << ','
      << c.cluster_radius_mm << ',' << c.cluster_radius_mm << ',' << c.component_keep_largest << ','
      << c.component_min_size << ',' << c.confirmed_disconnected_noise << ','
      << c.adaptive_local << ',' << c.local_knn << ','
      << c.local_pca_knn << ',' << c.local_mad_factor << ',' << c.local_regional_threshold << ','
      << c.local_end_threshold_multiplier << ',' << c.local_score_threshold << ','
      << output.injected_outlier_percent << ',' << output.injected_outlier_distance_mm << ','
      << output.injected_outlier_mode << ',' << d.points_before << ',' << d.points_after << ','
      << d.points_after_sor << ',' << d.points_after_ror << ',' << d.points_after_component << ','
      << d.points_after_local << ',' << d.removed_points << ',' << d.removed_percent << ','
      << d.removed_sor_points << ',' << d.removed_ror_points << ',' << d.removed_cluster_points << ','
      << d.removed_local_points << ','
      << d.top_before << ',' << d.top_after << ',' << d.top_removed_percent << ','
      << d.middle_before << ',' << d.middle_after << ',' << d.middle_removed_percent << ','
      << d.bottom_before << ',' << d.bottom_after << ',' << d.bottom_removed_percent << ','
      << d.length_before_mm << ',' << d.length_after_mm << ',' << d.length_change_mm << ','
      << d.length_change_percent << ',' << d.confirmed_component_length_change_mm << ','
      << d.unexplained_length_change_mm << ',' << d.possible_over_denoise << ','
      << d.warning_end_damage << ',' << d.invalid_over_filtering << ','
      << d.number_of_clusters << ',' << d.number_of_clusters << ','
      << d.largest_cluster_points << ',' << d.largest_cluster_points << ','
      << d.second_cluster_points << ',' << d.removed_cluster_points << ','
      << d.local_bottom.density.median << ',' << d.local_bottom.density.mad << ','
      << d.local_bottom.density.threshold << ','
      << d.local_middle.density.median << ',' << d.local_middle.density.mad << ','
      << d.local_middle.density.threshold << ','
      << d.local_top.density.median << ',' << d.local_top.density.mad << ','
      << d.local_top.density.threshold << ','
      << r.ply_load_ms << ',' << r.downsample_ms << ',' << d.sor_ms << ',' << d.ror_ms << ','
      << d.clustering_ms << ',' << d.clustering_ms << ',' << d.local_ms << ','
      << d.denoise_total_ms << ',' << d.denoise_total_ms << ','
      << r.alpha_mm << ',' << r.offset_mm << ','
      << r.alpha_wrap_ms << ',' << r.component_cleanup_ms << ',' << r.validation_ms << ',' << r.volume_ms << ','
      << r.online_core_ms << ',' << r.total_ms << ',' << r.total_ms << ','
      << v.vertices << ',' << v.faces << ','
      << v.connected_components << ',' << v.is_closed << ',' << v.is_closed << ',' << v.manifold << ','
      << v.self_intersection << ',' << v.outward_oriented << ',' << v.two_sided_wrap << ','
      << r.volume_mm3 << ','
      << r.volume_mm3 / 1000.0 << ',' << v.surface_area_mm2 << ',' << r.valid << ','
      << clean_csv(r.failure_reason) << '\n';
}

void write_density_csv(const fs::path& path, const RunOutput& output) {
  const DensityAnalysis density = analyze_nearest_neighbor_density(output.denoised_points,
      output.denoise_config.metric_end_fraction);
  fs::create_directories(path.parent_path().empty() ? fs::path(".") : path.parent_path());
  const bool needs_header = !fs::exists(path) || fs::file_size(path) == 0;
  std::ofstream out(path, std::ios::app);
  if (needs_header)
    out << "voxel_size_mm,denoise_method,end_protection,region,points,mean,median,std,p5,p25,p75,p95,density_runtime_ms\n";
  const auto row = [&](const char* region, std::size_t count, const DistributionStats& s) {
    out << std::setprecision(12) << output.result.voxel_size_mm << ','
        << denoise_method_name(output.denoise_config.method) << ','
        << output.denoise_config.end_protection_fraction << ',' << region << ',' << count << ','
        << s.mean << ',' << s.median << ',' << s.std << ',' << s.p5 << ',' << s.p25 << ','
        << s.p75 << ',' << s.p95 << ',' << density.runtime_ms << '\n';
  };
  row("OVERALL", density.overall_points, density.overall);
  row("BOTTOM", density.bottom_points, density.bottom);
  row("MIDDLE", density.middle_points, density.middle);
  row("TOP", density.top_points, density.top);
}

std::string determine_failure(const BenchmarkResult& r) {
  std::vector<std::string> reasons;
  const auto add = [&](const bool condition, const std::string& text) {
    if (condition) reasons.push_back(text);
  };
  add(r.validation.faces == 0, "EMPTY_MESH");
  add(!r.validation.is_closed, "NOT_WATERTIGHT");
  add(!r.validation.is_triangle_mesh, "NOT_TRIANGLE_MESH");
  add(r.validation.connected_components != 1, "MULTIPLE_COMPONENTS");
  add(!r.validation.manifold, "NON_MANIFOLD");
  add(r.validation.self_intersection, "SELF_INTERSECTION");
  add(!r.validation.outward_oriented, "NOT_OUTWARD_ORIENTED");
  add(r.validation.two_sided_wrap, "INVALID_TWO_SIDED_WRAP");
  add(!std::isfinite(r.volume_mm3) || r.volume_mm3 <= 0, "INVALID_VOLUME");
  std::ostringstream out;
  for (std::size_t i = 0; i < reasons.size(); ++i) {
    if (i) out << '|';
    out << reasons[i];
  }
  return out.str();
}

std::vector<Point> inject_outliers(const std::vector<Point>& input, const DatasetStats& stats,
                                   const double percent, const double distance_mm,
                                   const std::string& mode, const double voxel_mm) {
  if (!(percent > 0)) return input;
  std::vector<Point> output = input;
  const std::size_t count = std::max<std::size_t>(1, static_cast<std::size_t>(
      std::llround(percent * 0.01 * static_cast<double>(input.size()))));
  output.reserve(input.size() + count);
  const double z_span = stats.max[2] - stats.min[2];
  for (std::size_t i = 0; i < count; ++i) {
    const bool clustered = mode == "clustered" || mode == "hybrid";
    const std::size_t seed = clustered ? i / 8 : i;
    const std::size_t member = clustered ? i % 8 : 0;
    const double t1 = 0.1 + 0.8 * static_cast<double>((seed * 2654435761ULL) % 1000003ULL) / 1000003.0;
    const double t2 = 0.1 + 0.8 * static_cast<double>((seed * 2246822519ULL + 3266489917ULL) % 1000033ULL) / 1000033.0;
    const double z = stats.min[2] + t2 * z_span;
    const double jitter = clustered ? std::max(0.1, 0.22 * voxel_mm) : 0.0;
    const double angle = 0.7853981633974483 * static_cast<double>(member);
    const double j1 = jitter * std::cos(angle);
    const double j2 = jitter * std::sin(angle);
    switch (seed % 4) {
      case 0: output.emplace_back(stats.min[0] - distance_mm + j1, stats.min[1] + t1 * (stats.max[1] - stats.min[1]) + j2, z + j2); break;
      case 1: output.emplace_back(stats.max[0] + distance_mm + j1, stats.min[1] + t1 * (stats.max[1] - stats.min[1]) + j2, z + j2); break;
      case 2: output.emplace_back(stats.min[0] + t1 * (stats.max[0] - stats.min[0]) + j2, stats.min[1] - distance_mm + j1, z + j2); break;
      default: output.emplace_back(stats.min[0] + t1 * (stats.max[0] - stats.min[0]) + j2, stats.max[1] + distance_mm + j1, z + j2); break;
    }
  }
  return output;
}

RunOutput run_once(const std::vector<Point>& raw, const DatasetStats& raw_stats,
                   const double load_ms, const double voxel_mm,
                   const std::size_t max_points, const double alpha_mm,
                   const double offset_mm, const int run_index,
                   const DenoiseConfig& denoise_config = {},
                   const double outlier_percent = 0,
                   const double outlier_distance_mm = 0,
                   const std::string& outlier_mode = "dispersed",
                   const bool filter_only = false) {
  RunOutput output;
  output.denoise_config = denoise_config;
  output.injected_outlier_percent = outlier_percent;
  output.injected_outlier_distance_mm = outlier_distance_mm;
  output.injected_outlier_mode = outlier_mode;
  BenchmarkResult& r = output.result;
  r.voxel_size_mm = voxel_mm;
  r.original_points = raw.size();
  r.alpha_mm = alpha_mm;
  r.offset_mm = offset_mm;
  r.run_index = run_index;
  r.ply_load_ms = load_ms;

  std::vector<Point> processed;
  const std::vector<Point>* input = &raw;
  const auto down_begin = Clock::now();
  if (voxel_mm > 0) {
    processed = voxel_downsample_centroid(raw, voxel_mm, raw_stats);
    input = &processed;
  }
  if (max_points > 0 && max_points < input->size()) {
    processed = systematic_sample(*input, max_points);
    input = &processed;
  }
  const auto down_end = Clock::now();
  r.downsample_ms = elapsed_ms(down_begin, down_end);
  if (denoise_config.capture_stages) output.voxel_points = *input;
  std::vector<Point> protection_reference;
  std::vector<Point> injected;
  if (outlier_percent > 0) {
    const DatasetStats pre_denoise_stats = compute_dataset_stats(*input);
    protection_reference = *input;
    injected = inject_outliers(*input, pre_denoise_stats, outlier_percent,
                               outlier_distance_mm, outlier_mode, voxel_mm);
    input = &injected;
  }
  output.points_before_denoise = input->size();
  if (denoise_config.method != DenoiseMethod::none || denoise_config.clustering ||
      denoise_config.adaptive_local) {
    DenoiseOutput denoised = denoise_points(*input, denoise_config,
        protection_reference.empty() ? nullptr : &protection_reference);
    output.denoise_metrics = denoised.metrics;
    output.removed_points = std::move(denoised.removed);
    output.denoised_points = std::move(denoised.points);
    output.stages = std::move(denoised.stages);
    output.components = std::move(denoised.components);
    output.component_labels = std::move(denoised.component_labels);
    input = &output.denoised_points;
    r.denoise_ms = output.denoise_metrics.denoise_total_ms;
  } else {
    DenoiseConfig metric_config = denoise_config;
    DenoiseOutput metrics_only = denoise_points(*input, metric_config,
        protection_reference.empty() ? nullptr : &protection_reference);
    output.denoise_metrics = metrics_only.metrics;
    output.denoise_metrics.denoise_total_ms = 0;
    output.denoise_metrics.sor_ms = output.denoise_metrics.ror_ms =
        output.denoise_metrics.clustering_ms = output.denoise_metrics.local_ms = 0;
    output.denoised_points = std::move(metrics_only.points);
    input = &output.denoised_points;
  }
  r.input_points = input->size();

  if (filter_only) {
    if (output.denoise_metrics.invalid_over_filtering)
      r.failure_reason = "INVALID_OVER_FILTERING";
    r.valid = r.failure_reason.empty();
    r.reconstruction_only_ms = 0;
    r.online_core_ms = r.downsample_ms + r.denoise_ms;
    r.total_ms = r.online_core_ms;
    return output;
  }

  try {
    const auto wrap_begin = Clock::now();
    output.mesh = run_alpha_wrap(*input, alpha_mm, offset_mm);
    const auto wrap_end = Clock::now();
    r.alpha_wrap_ms = elapsed_ms(wrap_begin, wrap_end);

    const auto cleanup_begin = Clock::now();
    r.cleanup = cleanup_tiny_components(output.mesh);
    const auto cleanup_end = Clock::now();
    r.component_cleanup_ms = elapsed_ms(cleanup_begin, cleanup_end);

    const auto validation_begin = Clock::now();
    const DatasetStats input_stats = compute_dataset_stats(*input);
    r.validation = validate_mesh(output.mesh, *input, input_stats, offset_mm);
    const auto validation_end = Clock::now();
    r.validation_ms = elapsed_ms(validation_begin, validation_end);

    const auto volume_begin = Clock::now();
    const auto [signed_volume, absolute_volume] = mesh_volume(output.mesh);
    const auto volume_end = Clock::now();
    r.signed_volume_mm3 = signed_volume;
    r.volume_mm3 = absolute_volume;
    r.volume_ms = elapsed_ms(volume_begin, volume_end);
    r.failure_reason = determine_failure(r);
    if (output.denoise_metrics.invalid_over_filtering) {
      if (!r.failure_reason.empty()) r.failure_reason += '|';
      r.failure_reason += "INVALID_OVER_FILTERING";
    }
    r.valid = r.failure_reason.empty();
  } catch (const std::exception& e) {
    r.failure_reason = std::string("EXCEPTION:") + e.what();
  }
  r.reconstruction_only_ms = r.alpha_wrap_ms + r.volume_ms;
  r.online_core_ms = r.downsample_ms + r.denoise_ms + r.alpha_wrap_ms +
                     r.component_cleanup_ms + r.volume_ms;
  r.total_ms = r.online_core_ms + r.validation_ms;
  return output;
}

void print_result(const BenchmarkResult& r) {
  const auto& v = r.validation;
  std::cout << std::fixed << std::setprecision(3)
            << "RESULT run=" << r.run_index << " voxel=" << r.voxel_size_mm
            << " points=" << r.input_points << " alpha=" << r.alpha_mm
            << " offset=" << r.offset_mm << " downsample_ms=" << r.downsample_ms
            << " denoise_ms=" << r.denoise_ms
            << " alpha_wrap_ms=" << r.alpha_wrap_ms << " cleanup_ms=" << r.component_cleanup_ms
            << " validation_ms=" << r.validation_ms
            << " volume_ms=" << r.volume_ms << " online_core_ms=" << r.online_core_ms
            << " volume_ml=" << r.volume_mm3 / 1000.0 << " faces=" << v.faces
            << " closed=" << v.is_closed << " manifold=" << v.manifold
            << " raw_components=" << r.cleanup.raw_connected_components
            << " components=" << v.connected_components << " two_sided=" << v.two_sided_wrap
            << " valid=" << r.valid << " failure=" << r.failure_reason << "\n";
}

void write_meshes(const std::string& prefix, const Mesh& mesh) {
  if (prefix.empty() || mesh.is_empty()) return;
  const fs::path base(prefix);
  if (!base.parent_path().empty()) fs::create_directories(base.parent_path());
  const fs::path ply = base.string() + ".ply";
  const fs::path stl = base.string() + ".stl";
  if (!CGAL::IO::write_polygon_mesh(ply.string(), mesh, CGAL::parameters::stream_precision(17)))
    throw std::runtime_error("Failed to write " + ply.string());
  if (!CGAL::IO::write_polygon_mesh(stl.string(), mesh, CGAL::parameters::stream_precision(17)))
    throw std::runtime_error("Failed to write " + stl.string());
}

void run_voxel_benchmark(const std::vector<Point>& raw, const DatasetStats& stats,
                         const fs::path& csv_path, const std::vector<double>& sizes) {
  fs::create_directories(csv_path.parent_path().empty() ? fs::path(".") : csv_path.parent_path());
  std::ofstream csv(csv_path);
  csv << "voxel_size_mm,original_points,downsampled_points,downsampling_time_ms\n";
  csv << "Original," << raw.size() << ',' << raw.size() << ",0\n";
  for (const double voxel : sizes) {
    const auto begin = Clock::now();
    const auto sampled = voxel_downsample_centroid(raw, voxel, stats);
    const auto end = Clock::now();
    const double ms = elapsed_ms(begin, end);
    csv << std::setprecision(12) << voxel << ',' << raw.size() << ',' << sampled.size() << ',' << ms << '\n';
    std::cout << "VOXEL voxel=" << voxel << " points=" << sampled.size() << " ms=" << ms << '\n';
  }
}

void print_help() {
  std::cout <<
      "cucumber_alpha_wrap --input file.ply [mode/options]\n"
      "  --info\n"
      "  --voxel-benchmark-csv out.csv [--voxel-list 0.5,0.75,...]\n"
      "  --voxel mm --alpha mm --offset mm [--max-points N]\n"
      "  [--denoise none|sor|ror|sor_ror --sor-k K --sor-std R]\n"
      "  [--ror-radius mm --ror-min N --end-protection fraction --clustering]\n"
      "  [--cluster-radius mm --component-strategy largest|min_size --component-min-size N]\n"
      "  [--adaptive-local --local-knn K --local-pca-knn K --local-mad-factor F]\n"
      "  [--local-global-threshold --local-end-multiplier F --local-score-threshold F]\n"
      "  [--cleaning-output-dir path --component-csv out.csv]\n"
      "  [--repeat N] [--result-csv out.csv] [--output-prefix path]\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments args = parse_arguments(argc, argv);
    if (args.has("--help")) { print_help(); return 0; }
    const fs::path input_path = args.has("--input") ? fs::path(args.get("--input")) : find_default_input();
    const auto load_begin = Clock::now();
    const std::vector<Point> raw = read_ply_xyz(input_path);
    const auto load_end = Clock::now();
    const double load_ms = elapsed_ms(load_begin, load_end);
    const DatasetStats raw_stats = compute_dataset_stats(raw);
    print_stats(raw_stats, load_ms);
    if (args.has("--info")) return 0;

    if (args.has("--voxel-benchmark-csv")) {
      const auto sizes = parse_doubles(args.get("--voxel-list", "0.5,0.75,1,1.25,1.5,1.75,2,2.5,3"));
      run_voxel_benchmark(raw, raw_stats, args.get("--voxel-benchmark-csv"), sizes);
      return 0;
    }

    const double voxel_mm = std::stod(args.get("--voxel", "2.0"));
    const double alpha_mm = std::stod(args.get("--alpha", "30.0"));
    const double offset_mm = std::stod(args.get("--offset", "0.5"));
    const std::size_t max_points = static_cast<std::size_t>(std::stoull(args.get("--max-points", "0")));
    DenoiseConfig denoise_config;
    denoise_config.method = parse_denoise_method(args.get("--denoise", "none"));
    denoise_config.sor_k = std::stoi(args.get("--sor-k", "30"));
    denoise_config.sor_std_ratio = std::stod(args.get("--sor-std", "2.0"));
    denoise_config.ror_radius_mm = std::stod(args.get("--ror-radius", std::to_string(std::max(1.0, 2.5 * voxel_mm))));
    denoise_config.ror_min_neighbors = std::stoi(args.get("--ror-min", "3"));
    denoise_config.end_protection_fraction = std::stod(args.get("--end-protection", "0"));
    denoise_config.metric_end_fraction = std::stod(args.get("--end-region-fraction",
        denoise_config.end_protection_fraction > 0 ? args.get("--end-protection") : "0.05"));
    denoise_config.clustering = args.has("--clustering");
    denoise_config.cluster_radius_mm = std::stod(args.get("--cluster-radius", std::to_string(std::max(1.0, 3.0 * voxel_mm))));
    denoise_config.component_keep_largest = args.get("--component-strategy", "largest") != "min_size";
    denoise_config.component_min_size = static_cast<std::size_t>(
        std::stoull(args.get("--component-min-size", "1")));
    denoise_config.component_end_keep_distance_multiplier = std::stod(
        args.get("--component-end-distance-multiplier", "2.0"));
    denoise_config.confirmed_disconnected_noise =
        args.has("--confirmed-disconnected-noise");
    denoise_config.adaptive_local = args.has("--adaptive-local");
    denoise_config.local_knn = std::stoi(args.get("--local-knn", "20"));
    denoise_config.local_pca_knn = std::stoi(args.get("--local-pca-knn", "20"));
    denoise_config.local_mad_factor = std::stod(args.get("--local-mad-factor", "5.0"));
    denoise_config.local_regional_threshold = !args.has("--local-global-threshold");
    denoise_config.local_end_threshold_multiplier = std::stod(
        args.get("--local-end-multiplier", "1.5"));
    denoise_config.local_score_threshold = std::stod(
        args.get("--local-score-threshold", "0.8"));
    denoise_config.capture_stages = args.has("--cleaning-output-dir");
    const double outlier_percent = std::stod(args.get("--outlier-percent", "0"));
    const double outlier_distance = std::stod(args.get("--outlier-distance", "0"));
    const std::string outlier_mode = args.get("--outlier-mode", "dispersed");
    if (outlier_mode != "dispersed" && outlier_mode != "clustered" && outlier_mode != "hybrid")
      throw std::runtime_error("--outlier-mode must be dispersed, clustered, or hybrid");
    const int repeats = std::max(1, std::stoi(args.get("--repeat", "1")));
    const bool do_warmup = args.has("--repeat");
    const bool filter_only = args.has("--filter-only");
    const std::string csv_path = args.get("--result-csv");
    const std::string denoise_csv_path = args.get("--denoise-result-csv");
    RunOutput last;
    if (do_warmup) {
      std::cout << "Warm-up (excluded from CSV/statistics)\n";
      const RunOutput warmup = run_once(raw, raw_stats, load_ms, voxel_mm, max_points,
                                        alpha_mm, offset_mm, -1, denoise_config,
                                        outlier_percent, outlier_distance, outlier_mode,
                                        filter_only);
      print_result(warmup.result);
    }
    for (int i = 0; i < repeats; ++i) {
      last = run_once(raw, raw_stats, load_ms, voxel_mm, max_points,
                      alpha_mm, offset_mm, i, denoise_config,
                      outlier_percent, outlier_distance, outlier_mode, filter_only);
      print_result(last.result);
      if (!csv_path.empty()) write_result_csv(csv_path, last.result);
      if (!denoise_csv_path.empty()) write_denoise_result_csv(denoise_csv_path, last);
    }
    write_meshes(args.get("--output-prefix"), last.mesh);
    const std::string point_prefix = args.get("--output-points-prefix");
    if (!point_prefix.empty()) {
      write_xyz_ply(point_prefix + "_remaining.ply", last.denoised_points);
      write_xyz_ply(point_prefix + "_removed.ply", last.removed_points);
    }
    const std::string cleaning_dir = args.get("--cleaning-output-dir");
    if (!cleaning_dir.empty()) {
      const fs::path dir(cleaning_dir);
      fs::create_directories(dir);
      write_xyz_ply(dir / "00_original.ply", raw);
      write_xyz_ply(dir / "01_voxel.ply", last.voxel_points);
      write_xyz_ply(dir / "02_sor_remaining.ply", last.stages.after_sor);
      write_xyz_ply(dir / "02_sor_removed.ply", last.stages.removed_sor);
      write_xyz_ply(dir / "03_ror_remaining.ply", last.stages.after_ror);
      write_xyz_ply(dir / "03_ror_removed.ply", last.stages.removed_ror);
      write_xyz_ply(dir / "04_component_remaining.ply", last.stages.after_component);
      write_xyz_ply(dir / "04_component_removed.ply", last.stages.removed_component);
      write_xyz_ply(dir / "removed_clusters.ply", last.stages.removed_component);
      write_xyz_ply(dir / "05_local_remaining.ply", last.stages.after_local);
      write_xyz_ply(dir / "05_local_removed.ply", last.stages.removed_local);
      write_xyz_ply(dir / "removed_local_outliers.ply", last.stages.removed_local);
      write_xyz_ply(dir / "06_final_cleaned.ply", last.denoised_points);
      if (!last.component_labels.empty())
        write_component_rgb_ply(dir / "component_colored.ply", last.stages.after_ror,
                                last.component_labels, last.components);
    }
    if (args.has("--component-csv"))
      write_component_csv(args.get("--component-csv"), last.components);
    if (args.has("--density-csv")) write_density_csv(args.get("--density-csv"), last);
    return last.result.valid ? 0 : 2;
  } catch (const std::exception& e) {
    std::cerr << "FATAL: " << e.what() << '\n';
    return 1;
  }
}
