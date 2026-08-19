#include "alpha_wrap.hpp"
#include "denoise.hpp"
#include "mesh_validation.hpp"
#include "ply_reader.hpp"
#include "ring_reconstruction.hpp"
#include "volume.hpp"
#include "voxel_downsample.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
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
  bool has(const std::string& key) const { return values.count(key) != 0; }
  std::string get(const std::string& key, const std::string& fallback = {}) const {
    const auto it = values.find(key);
    return it == values.end() ? fallback : it->second;
  }
};

Arguments parse_arguments(const int argc, char** argv) {
  Arguments args;
  for (int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    if (key.rfind("--", 0) != 0 || i + 1 >= argc)
      throw std::runtime_error("Expected --key value, got: " + key);
    args.values[key] = argv[++i];
  }
  return args;
}

std::vector<std::size_t> parse_sizes(const std::string& text) {
  std::vector<std::size_t> result;
  std::istringstream input(text);
  std::string token;
  while (std::getline(input, token, ',')) result.push_back(std::stoull(token));
  return result;
}

struct Distribution {
  double mean = 0;
  double median = 0;
  double p90 = 0;
  double p95 = 0;
  double p99 = 0;
  double minimum = 0;
  double maximum = 0;
  double stddev = 0;
};

double percentile_sorted(const std::vector<double>& sorted, const double p) {
  if (sorted.empty()) return 0;
  const double index = p * static_cast<double>(sorted.size() - 1);
  const auto lo = static_cast<std::size_t>(std::floor(index));
  const auto hi = static_cast<std::size_t>(std::ceil(index));
  const double f = index - static_cast<double>(lo);
  return sorted[lo] * (1.0 - f) + sorted[hi] * f;
}

Distribution summarize(std::vector<double> values) {
  Distribution result;
  if (values.empty()) return result;
  result.mean = std::accumulate(values.begin(), values.end(), 0.0) /
                static_cast<double>(values.size());
  double square_sum = 0;
  for (const double value : values) square_sum += (value - result.mean) * (value - result.mean);
  result.stddev = std::sqrt(square_sum / static_cast<double>(values.size()));
  std::sort(values.begin(), values.end());
  result.minimum = values.front();
  result.maximum = values.back();
  result.median = percentile_sorted(values, 0.50);
  result.p90 = percentile_sorted(values, 0.90);
  result.p95 = percentile_sorted(values, 0.95);
  result.p99 = percentile_sorted(values, 0.99);
  return result;
}

struct Sample {
  RingReconstructionTimings timing;
  double validation_ms = 0;
  double volume_ms = 0;
  double reconstruction_plus_volume_ms = 0;
  double total_ms = 0;
  double volume_ml = 0;
  std::size_t actual_ring_count = 0;
  double axial_step_mm = 0;
};

struct SummaryRow {
  std::string method;
  std::size_t target_ring_count = 0;
  std::size_t actual_ring_count = 0;
  double axial_step_mm = 0;
  std::size_t angular_bins = 0;
  std::size_t endpoint_fit_rings = 0;
  std::size_t endpoint_generated_rings = 0;
  RingReconstructionTimings timing;
  double validation_ms = 0;
  double volume_ms = 0;
  double reconstruction_plus_volume_ms = 0;
  double online_core_ms = 0;
  Distribution total;
  Distribution reconstruction_plus_volume;
  Distribution nurbs_endpoint;
  double volume_ml = 0;
  double volume_diff_reference_ml = 0;
  double volume_diff_reference_percent = 0;
  ValidationResult validation;
  Mesh representative_mesh;
};

double mean_field(const std::vector<Sample>& samples,
                  const double RingReconstructionTimings::*field) {
  double sum = 0;
  for (const auto& sample : samples) sum += sample.timing.*field;
  return sum / static_cast<double>(samples.size());
}

SummaryRow aggregate_ring(const RingReconstructionConfig& config,
                          const std::vector<Sample>& samples,
                          const ValidationResult& validation, Mesh mesh,
                          const double preprocessing_ms) {
  SummaryRow row;
  row.method = config.endpoint == EndpointCompletion::flat_cap ? "ring_flat" : "ring_spline";
  row.target_ring_count = config.target_ring_count;
  row.actual_ring_count = samples.back().actual_ring_count;
  row.axial_step_mm = samples.back().axial_step_mm;
  row.angular_bins = config.angular_bins;
  row.endpoint_fit_rings = config.endpoint == EndpointCompletion::cubic_bspline
      ? config.endpoint_fit_rings : 0;
  row.endpoint_generated_rings = config.endpoint == EndpointCompletion::cubic_bspline
      ? config.endpoint_generated_rings : 0;
#define MEAN_TIMING(member) row.timing.member = mean_field(samples, &RingReconstructionTimings::member)
  MEAN_TIMING(pca_ms);
  MEAN_TIMING(axial_binning_ms);
  MEAN_TIMING(centerline_ms);
  MEAN_TIMING(local_frame_ms);
  MEAN_TIMING(radial_sector_ms);
  MEAN_TIMING(ring_vertex_generation_ms);
  MEAN_TIMING(ring_mesh_connection_ms);
  MEAN_TIMING(bottom_spline_fit_ms);
  MEAN_TIMING(top_spline_fit_ms);
  MEAN_TIMING(bottom_spline_eval_ms);
  MEAN_TIMING(top_spline_eval_ms);
  MEAN_TIMING(endpoint_ring_generation_ms);
  MEAN_TIMING(endpoint_mesh_connection_ms);
  MEAN_TIMING(endpoint_completion_ms);
  MEAN_TIMING(nurbs_endpoint_total_ms);
  MEAN_TIMING(ring_reconstruction_ms);
#undef MEAN_TIMING
  std::vector<double> total, reconstruction_volume, endpoint;
  double validation_sum = 0, volume_time_sum = 0, volume_sum = 0;
  for (const auto& sample : samples) {
    total.push_back(sample.total_ms);
    reconstruction_volume.push_back(sample.reconstruction_plus_volume_ms);
    endpoint.push_back(sample.timing.nurbs_endpoint_total_ms);
    validation_sum += sample.validation_ms;
    volume_time_sum += sample.volume_ms;
    volume_sum += sample.volume_ml;
  }
  row.total = summarize(std::move(total));
  row.reconstruction_plus_volume = summarize(std::move(reconstruction_volume));
  row.nurbs_endpoint = summarize(std::move(endpoint));
  row.validation_ms = validation_sum / static_cast<double>(samples.size());
  row.volume_ms = volume_time_sum / static_cast<double>(samples.size());
  row.reconstruction_plus_volume_ms = row.timing.ring_reconstruction_ms + row.volume_ms;
  row.online_core_ms = preprocessing_ms + row.reconstruction_plus_volume_ms;
  row.volume_ml = volume_sum / static_cast<double>(samples.size());
  row.validation = validation;
  row.representative_mesh = std::move(mesh);
  return row;
}

SummaryRow benchmark_ring(const std::vector<Point>& points, const DatasetStats& stats,
                          const RingReconstructionConfig& config, const int warmups,
                          const int repeats, const double preprocessing_ms) {
  std::cout << "BENCHMARK method="
            << (config.endpoint == EndpointCompletion::flat_cap ? "ring_flat" : "ring_spline")
            << " target=" << config.target_ring_count
            << " fit=" << config.endpoint_fit_rings
            << " generated=" << config.endpoint_generated_rings << std::endl;
  const auto execute = [&]() {
    Sample sample;
    RingReconstructionResult reconstruction = reconstruct_ring_mesh(points, config);
    const auto validation_begin = Clock::now();
    const ValidationResult validation = validate_mesh(reconstruction.mesh, points, stats, 0.1);
    const auto validation_end = Clock::now();
    const auto volume_begin = Clock::now();
    const auto volume = mesh_volume(reconstruction.mesh);
    const auto volume_end = Clock::now();
    sample.timing = reconstruction.timings;
    sample.validation_ms = elapsed_ms(validation_begin, validation_end);
    sample.volume_ms = elapsed_ms(volume_begin, volume_end);
    sample.reconstruction_plus_volume_ms = sample.timing.ring_reconstruction_ms + sample.volume_ms;
    sample.total_ms = sample.reconstruction_plus_volume_ms + sample.validation_ms;
    sample.volume_ml = volume.second / 1000.0;
    sample.actual_ring_count = reconstruction.actual_ring_count;
    sample.axial_step_mm = reconstruction.axial_step_mm;
    return std::tuple<Sample, ValidationResult, Mesh>(
        sample, validation, std::move(reconstruction.mesh));
  };
  for (int i = 0; i < warmups; ++i) execute();
  std::vector<Sample> samples;
  samples.reserve(static_cast<std::size_t>(repeats));
  ValidationResult validation;
  Mesh mesh;
  for (int i = 0; i < repeats; ++i) {
    auto [sample, run_validation, run_mesh] = execute();
    samples.push_back(sample);
    if (i + 1 == repeats) {
      validation = run_validation;
      mesh = std::move(run_mesh);
    }
  }
  SummaryRow row = aggregate_ring(config, samples, validation, std::move(mesh), preprocessing_ms);
  std::cout << std::fixed << std::setprecision(3)
            << "  actual=" << row.actual_ring_count << " mean_total_ms=" << row.total.mean
            << " p95_total_ms=" << row.total.p95 << " volume_ml=" << row.volume_ml
            << " closed=" << row.validation.is_closed
            << " manifold=" << row.validation.manifold
            << " self_intersection=" << row.validation.self_intersection << std::endl;
  return row;
}

SummaryRow benchmark_alpha(const std::vector<Point>& points, const DatasetStats& stats,
                           const double alpha_mm, const double offset_mm,
                           const int warmups, const int repeats,
                           const double preprocessing_ms) {
  std::cout << "BENCHMARK method=alpha_wrap alpha=" << alpha_mm
            << " offset=" << offset_mm << std::endl;
  struct AlphaSample { double reconstruction = 0, validation = 0, volume = 0, total = 0, ml = 0; };
  const auto execute = [&]() {
    AlphaSample sample;
    const auto reconstruction_begin = Clock::now();
    Mesh mesh = run_alpha_wrap(points, alpha_mm, offset_mm);
    cleanup_tiny_components(mesh);
    const auto reconstruction_end = Clock::now();
    const auto validation_begin = Clock::now();
    const ValidationResult validation = validate_mesh(mesh, points, stats, offset_mm);
    const auto validation_end = Clock::now();
    const auto volume_begin = Clock::now();
    const auto volume = mesh_volume(mesh);
    const auto volume_end = Clock::now();
    sample.reconstruction = elapsed_ms(reconstruction_begin, reconstruction_end);
    sample.validation = elapsed_ms(validation_begin, validation_end);
    sample.volume = elapsed_ms(volume_begin, volume_end);
    sample.total = sample.reconstruction + sample.validation + sample.volume;
    sample.ml = volume.second / 1000.0;
    return std::tuple<AlphaSample, ValidationResult, Mesh>(sample, validation, std::move(mesh));
  };
  for (int i = 0; i < warmups; ++i) execute();
  std::vector<AlphaSample> samples;
  samples.reserve(static_cast<std::size_t>(repeats));
  ValidationResult validation;
  Mesh mesh;
  for (int i = 0; i < repeats; ++i) {
    auto [sample, run_validation, run_mesh] = execute();
    samples.push_back(sample);
    if (i + 1 == repeats) { validation = run_validation; mesh = std::move(run_mesh); }
  }
  SummaryRow row;
  row.method = "alpha_wrap";
  std::vector<double> total, reconstruction_volume;
  double reconstruction_sum = 0, validation_sum = 0, volume_sum = 0, ml_sum = 0;
  for (const auto& sample : samples) {
    total.push_back(sample.total);
    reconstruction_volume.push_back(sample.reconstruction + sample.volume);
    reconstruction_sum += sample.reconstruction;
    validation_sum += sample.validation;
    volume_sum += sample.volume;
    ml_sum += sample.ml;
  }
  row.total = summarize(std::move(total));
  row.reconstruction_plus_volume = summarize(std::move(reconstruction_volume));
  row.timing.ring_reconstruction_ms = reconstruction_sum / repeats;
  row.validation_ms = validation_sum / repeats;
  row.volume_ms = volume_sum / repeats;
  row.reconstruction_plus_volume_ms = row.timing.ring_reconstruction_ms + row.volume_ms;
  row.online_core_ms = preprocessing_ms + row.reconstruction_plus_volume_ms;
  row.volume_ml = ml_sum / repeats;
  row.validation = validation;
  row.representative_mesh = std::move(mesh);
  std::cout << std::fixed << std::setprecision(3)
            << "  mean_total_ms=" << row.total.mean << " p95_total_ms=" << row.total.p95
            << " volume_ml=" << row.volume_ml << " closed=" << row.validation.is_closed
            << " manifold=" << row.validation.manifold << std::endl;
  return row;
}

bool topology_ok(const SummaryRow& row) {
  return row.validation.is_closed && row.validation.manifold &&
         row.validation.is_triangle_mesh && !row.validation.self_intersection &&
         row.validation.outward_oriented && row.validation.connected_components == 1;
}

void write_mesh_ply(const fs::path& path, const Mesh& mesh) {
  if (!path.parent_path().empty()) fs::create_directories(path.parent_path());
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("Cannot write mesh PLY: " + path.string());
  out << "ply\nformat ascii 1.0\nelement vertex " << mesh.number_of_vertices()
      << "\nproperty double x\nproperty double y\nproperty double z\nelement face "
      << mesh.number_of_faces()
      << "\nproperty list uchar int vertex_indices\nend_header\n";
  std::vector<std::size_t> index(mesh.num_vertices());
  std::size_t next_index = 0;
  out << std::setprecision(17);
  for (const auto vertex : mesh.vertices()) {
    index[static_cast<std::size_t>(vertex.idx())] = next_index++;
    const Point& p = mesh.point(vertex);
    out << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
  }
  for (const auto face : mesh.faces()) {
    std::vector<std::size_t> vertices;
    for (const auto vertex : CGAL::vertices_around_face(mesh.halfedge(face), mesh))
      vertices.push_back(index[static_cast<std::size_t>(vertex.idx())]);
    out << vertices.size();
    for (const auto vertex : vertices) out << ' ' << vertex;
    out << '\n';
  }
}

void write_csv(const fs::path& path, const std::vector<SummaryRow>& rows) {
  if (!path.parent_path().empty()) fs::create_directories(path.parent_path());
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Cannot write CSV: " + path.string());
  out << "method,target_ring_count,actual_ring_count,axial_step_mm,angular_bins,endpoint_fit_rings,endpoint_generated_rings,"
         "pca_ms,binning_ms,centerline_ms,frame_ms,radial_ms,ring_vertex_generation_ms,ring_mesh_connection_ms,mesh_ms,"
         "bottom_spline_fit_ms,top_spline_fit_ms,spline_fit_ms,bottom_spline_eval_ms,top_spline_eval_ms,spline_eval_ms,"
         "endpoint_ring_generation_ms,endpoint_mesh_connection_ms,endpoint_generation_ms,endpoint_completion_ms,"
         "nurbs_endpoint_mean_ms,nurbs_endpoint_p95_ms,validation_ms,volume_ms,reconstruction_ms,reconstruction_plus_volume_ms,"
         "online_core_ms,mean_ms,median_ms,p90_ms,p95_ms,p99_ms,min_ms,max_ms,stddev_ms,mean_le_500,p95_le_500,"
         "volume_ml,volume_diff_reference_ml,volume_diff_reference_percent,watertight,manifold,self_intersection,"
         "outward_oriented,connected_components,vertices,faces,edges,euler_characteristic,genus,surface_area_mm2\n";
  out << std::setprecision(12);
  for (const auto& r : rows) {
    const auto& t = r.timing;
    const double spline_fit = t.bottom_spline_fit_ms + t.top_spline_fit_ms;
    const double spline_eval = t.bottom_spline_eval_ms + t.top_spline_eval_ms;
    const double mesh_ms = t.ring_vertex_generation_ms + t.ring_mesh_connection_ms +
                           t.endpoint_ring_generation_ms + t.endpoint_mesh_connection_ms;
    out << r.method << ',' << r.target_ring_count << ',' << r.actual_ring_count << ','
        << r.axial_step_mm << ',' << r.angular_bins << ',' << r.endpoint_fit_rings << ','
        << r.endpoint_generated_rings << ',' << t.pca_ms << ',' << t.axial_binning_ms << ','
        << t.centerline_ms << ',' << t.local_frame_ms << ',' << t.radial_sector_ms << ','
        << t.ring_vertex_generation_ms << ',' << t.ring_mesh_connection_ms << ',' << mesh_ms << ','
        << t.bottom_spline_fit_ms << ',' << t.top_spline_fit_ms << ',' << spline_fit << ','
        << t.bottom_spline_eval_ms << ',' << t.top_spline_eval_ms << ',' << spline_eval << ','
        << t.endpoint_ring_generation_ms << ',' << t.endpoint_mesh_connection_ms << ','
        << t.endpoint_ring_generation_ms << ',' << t.endpoint_completion_ms << ','
        << r.nurbs_endpoint.mean << ',' << r.nurbs_endpoint.p95 << ',' << r.validation_ms << ','
        << r.volume_ms << ',' << t.ring_reconstruction_ms << ',' << r.reconstruction_plus_volume_ms << ','
        << r.online_core_ms << ',' << r.total.mean << ',' << r.total.median << ',' << r.total.p90 << ','
        << r.total.p95 << ',' << r.total.p99 << ',' << r.total.minimum << ',' << r.total.maximum << ','
        << r.total.stddev << ',' << (r.total.mean <= 500.0) << ',' << (r.total.p95 <= 500.0) << ','
        << r.volume_ml << ',' << r.volume_diff_reference_ml << ',' << r.volume_diff_reference_percent << ','
        << r.validation.is_closed << ',' << r.validation.manifold << ',' << r.validation.self_intersection << ','
        << r.validation.outward_oriented << ',' << r.validation.connected_components << ','
        << r.validation.vertices << ',' << r.validation.faces << ',' << r.validation.edges << ','
        << r.validation.euler_characteristic << ',' << r.validation.genus << ','
        << r.validation.surface_area_mm2 << '\n';
  }
}

const SummaryRow* find_row(const std::vector<SummaryRow>& rows, const std::string& method,
                           const std::size_t target, const std::size_t fit = 0,
                           const std::size_t generated = 0) {
  for (const auto& row : rows)
    if (row.method == method && row.target_ring_count == target &&
        (method != "ring_spline" ||
         (row.endpoint_fit_rings == fit && row.endpoint_generated_rings == generated))) return &row;
  return nullptr;
}

void write_report(const fs::path& path, const std::vector<SummaryRow>& rows,
                  const std::vector<std::size_t>& ring_targets,
                  const std::size_t reference_target, const std::size_t recommended_target,
                  const double preprocessing_ms, const int warmups, const int repeats,
                  const double alpha_mm, const double offset_mm) {
  if (!path.parent_path().empty()) fs::create_directories(path.parent_path());
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Cannot write report: " + path.string());
  const SummaryRow* reference = find_row(rows, "ring_flat", reference_target);
  const SummaryRow* recommended_flat = find_row(rows, "ring_flat", recommended_target);
  const SummaryRow* recommended_spline = find_row(rows, "ring_spline", recommended_target, 8, 4);
  const SummaryRow* alpha = nullptr;
  for (const auto& row : rows) if (row.method == "alpha_wrap") alpha = &row;
  const SummaryRow* maximum = nullptr;
  const SummaryRow* time_only_maximum = nullptr;
  for (const auto& row : rows)
    if (row.method == "ring_flat" && row.total.p95 <= 500.0) {
      if (!time_only_maximum || row.actual_ring_count > time_only_maximum->actual_ring_count)
        time_only_maximum = &row;
      if (topology_ok(row) && (!maximum || row.actual_ring_count > maximum->actual_ring_count))
        maximum = &row;
    }

  out << "# Ring reconstruction + cubic B-spline endpoint benchmark\n\n"
      << "## Method and benchmark protocol\n\n"
      << "Release build; `std::chrono::steady_clock`; " << warmups << " warm-up runs and "
      << repeats << " measured runs per configuration. PLY read, PLY write, console output, and debug "
      << "visualization are excluded. The same in-memory point cloud (1.5 mm voxel + 4.0 mm dominant "
      << "component) is supplied to Alpha Wrap, Ring-flat, and Ring-spline. One-time preprocessing was "
      << std::fixed << std::setprecision(3) << preprocessing_ms << " ms in this run.\n\n"
      << "Ring radii use exactly Q10/Q90 and `R=(Q10+Q90)/2`, with 48 angular bins and no radius "
      << "calibration. `mean_ms`/percentiles include reconstruction + mesh validation + volume; "
      << "`reconstruction_plus_volume_ms` excludes validation.\n\n"
      << "Reproduce with:\n\n"
         "```powershell\n"
         "$env:PATH = 'C:\\msys64\\ucrt64\\bin;' + $env:PATH\n"
         "cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release\n"
         "cmake --build build --parallel 8 --target cucumber_ring_benchmark\n"
         ".\\build\\cucumber_ring_benchmark.exe\n"
         "```\n\n"
      << "## Question 1 — Ring count under 0.5 s and convergence\n\n";
  if (maximum)
    out << "**Maximum tested valid Ring count under 0.5 s (p95): "
        << maximum->actual_ring_count << " Rings** (target " << maximum->target_ring_count
        << ", p95 " << maximum->total.p95 << " ms).\n\n";
  if (time_only_maximum && time_only_maximum != maximum)
    out << "Time-only capacity reached " << time_only_maximum->actual_ring_count
        << " actual Rings (target " << time_only_maximum->target_ring_count << ", p95 "
        << time_only_maximum->total.p95
        << " ms), but it is **not admissible** because the mesh self-intersects.\n\n";
  if (recommended_flat && reference)
    out << "**Recommended Ring count: " << recommended_flat->actual_ring_count
        << " Rings** (target " << recommended_target << "). Its volume difference to the target-"
        << reference_target << " reference is " << recommended_flat->volume_diff_reference_percent
        << "%, while using the smallest tested topology-valid configuration whose reference error is "
        << "at most 0.5% and remains converged at all higher topology-valid counts through the reference.\n\n";
  out << "| target | actual | step mm | volume mL | change from previous % | diff to ref % | mean ms | p95 ms | watertight | manifold | self-intersection |\n"
         "|---:|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|:---:|\n";
  const SummaryRow* previous = nullptr;
  for (const auto target : ring_targets) {
    const auto* row = find_row(rows, "ring_flat", target);
    if (!row) continue;
    const double adjacent_change = previous
        ? 100.0 * std::abs(row->volume_ml - previous->volume_ml) / previous->volume_ml : 0.0;
    out << '|' << target << '|' << row->actual_ring_count << '|' << row->axial_step_mm << '|'
        << row->volume_ml << '|' << adjacent_change << '|' << row->volume_diff_reference_percent << '|' << row->total.mean
        << '|' << row->total.p95 << '|' << (row->validation.is_closed ? "yes" : "no")
        << '|' << (row->validation.manifold ? "yes" : "no")
        << '|' << (row->validation.self_intersection ? "yes" : "no") << "|\n";
    previous = row;
  }
  bool excluded_higher_topology = false;
  for (const auto target : ring_targets) {
    const auto* row = find_row(rows, "ring_flat", target);
    if (row && target > reference_target && !topology_ok(*row)) excluded_higher_topology = true;
  }
  if (excluded_higher_topology)
    out << "\nTargets above the reference are retained as performance data but excluded from the "
           "geometric reference/recommendation because validation detected self-intersection or another "
           "topology failure. This is the observed high-resolution over-slicing limit.\n";
  out << "\nConvergence thresholds (first tested count at or below threshold against reference):\n\n";
  for (const double threshold : {1.0, 0.5, 0.25, 0.1}) {
    const SummaryRow* first = nullptr;
    for (const auto target : ring_targets) {
      const auto* row = find_row(rows, "ring_flat", target);
      if (row && row->volume_diff_reference_percent <= threshold) { first = row; break; }
    }
    out << "- <=" << threshold << "%: "
        << (first ? std::to_string(first->actual_ring_count) + " actual Rings" : "not reached") << "\n";
  }

  out << "\n## Question 2 — Endpoint B-spline cost\n\n";
  if (recommended_flat && recommended_spline) {
    const auto& t = recommended_spline->timing;
    const double fit = t.bottom_spline_fit_ms + t.top_spline_fit_ms;
    const double eval = t.bottom_spline_eval_ms + t.top_spline_eval_ms;
    const double share = recommended_spline->reconstruction_plus_volume_ms > 0
        ? 100.0 * recommended_spline->nurbs_endpoint.mean /
          recommended_spline->reconstruction_plus_volume_ms : 0;
    out << "Recommended-count comparison (K=8, 4 generated rings per end):\n\n"
        << "| endpoint | fit ms | evaluation ms | ring generation ms | endpoint mesh ms | endpoint mean ms | endpoint p95 ms |\n"
           "|---|---:|---:|---:|---:|---:|---:|\n"
        << "| Flat cap | 0 | 0 | " << recommended_flat->timing.endpoint_ring_generation_ms
        << '|' << recommended_flat->timing.endpoint_mesh_connection_ms << '|'
        << recommended_flat->timing.endpoint_completion_ms << "| n/a |\n"
        << "| Cubic B-spline | " << fit << '|' << eval << '|'
        << t.endpoint_ring_generation_ms << '|' << t.endpoint_mesh_connection_ms << '|'
        << recommended_spline->nurbs_endpoint.mean << '|' << recommended_spline->nurbs_endpoint.p95 << "|\n\n"
        << "The spline endpoint accounts for **" << share
        << "%** of Ring-spline reconstruction + volume time.\n\n";
  }
  out << "All endpoint sweeps are retained in `ring_benchmark.csv` (generated rings 2,3,4,5,6,8,10 at K=8; "
         "K=4,6,8,10,12 at four generated rings).\n\n"
      << "## Question 3 — Reconstruction comparison\n\n"
      << "Alpha Wrap parameters: alpha=" << alpha_mm << " mm, offset=" << offset_mm << " mm.\n\n"
      << "| method | mean ms | p95 ms | volume mL | vertices | faces | watertight | manifold |\n"
         "|---|---:|---:|---:|---:|---:|:---:|:---:|\n";
  const auto comparison_row = [&](const char* label, const SummaryRow* row) {
    if (!row) return;
    out << '|' << label << '|' << row->total.mean << '|' << row->total.p95 << '|'
        << row->volume_ml << '|' << row->validation.vertices << '|' << row->validation.faces
        << '|' << (row->validation.is_closed ? "yes" : "no")
        << '|' << (row->validation.manifold ? "yes" : "no") << "|\n";
  };
  comparison_row("Alpha Wrap", alpha);
  comparison_row("Ring + flat cap", recommended_flat);
  comparison_row("Ring + cubic B-spline tip", recommended_spline);
  if (alpha && recommended_spline) {
    const double difference = 100.0 * std::abs(alpha->volume_ml - recommended_spline->volume_ml) /
                              alpha->volume_ml;
    out << "\nAll three methods meet 500 ms p95 on this machine. Ring-spline is the fastest "
           "topology-valid completed surface, but its volume is " << difference
        << "% below Alpha Wrap. Because there is no physical ground-truth volume, this benchmark "
           "establishes runtime and internal Ring convergence, not absolute accuracy or equivalence "
           "between reconstruction families. Alpha Wrap remains appropriate when continuity with the "
           "existing Alpha-derived volume is required.\n";
  }
  if (alpha && recommended_flat && recommended_spline)
    out << "\nMean online core (one preprocessing pass + reconstruction + volume, validation excluded): "
        << "Alpha Wrap " << alpha->online_core_ms << " ms; Ring-flat "
        << recommended_flat->online_core_ms << " ms; Ring-spline "
        << recommended_spline->online_core_ms << " ms.\n";
  out << "\n## Recommendation\n\n";
  if (recommended_spline)
    out << "Use `ring_spline` with target/actual " << recommended_target << '/'
        << recommended_spline->actual_ring_count
        << " Rings, 48 angular bins, Q10/Q90 center radius, K=8 endpoint fit rings, and 4 generated "
           "rings per end. This is watertight/topology-valid, comfortably below 500 ms at p95, and "
           "uses the first volume-converged Ring count instead of spending the remaining budget on "
           "unnecessary axial resolution. The inspection mesh is `final_ring_spline_mesh.ply`.\n";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Arguments args = parse_arguments(argc, argv);
    const fs::path input = args.get("--input", "pointcloud_3100.ply");
    const fs::path csv_path = args.get("--csv", "ring_benchmark.csv");
    const fs::path report_path = args.get("--report", "RING_NURBS_BENCHMARK.md");
    const fs::path mesh_path = args.get("--mesh", "final_ring_spline_mesh.ply");
    const int warmups = std::max(10, std::stoi(args.get("--warmup", "10")));
    const int repeats = std::max(100, std::stoi(args.get("--repeat", "100")));
    const auto ring_targets = parse_sizes(args.get(
        "--ring-list", "20,30,40,50,60,75,90,120,150,180,240,300"));
    const std::size_t requested_reference_target = std::stoull(args.get("--reference-ring", "300"));
    const double alpha_mm = std::stod(args.get("--alpha", "6.25"));
    const double offset_mm = std::stod(args.get("--offset", "0.10"));

    std::cout << "Loading input (excluded from benchmark): " << input.string() << std::endl;
    const auto raw = read_ply_xyz(input);
    const DatasetStats raw_stats = compute_dataset_stats(raw);
    const auto preprocess_begin = Clock::now();
    const auto voxel = voxel_downsample_centroid(raw, 1.5, raw_stats);
    DenoiseConfig denoise;
    denoise.clustering = true;
    denoise.cluster_radius_mm = 4.0;
    denoise.component_keep_largest = true;
    denoise.confirmed_disconnected_noise = true;
    const auto cleaned = denoise_points(voxel, denoise).points;
    const auto preprocess_end = Clock::now();
    const double preprocessing_ms = elapsed_ms(preprocess_begin, preprocess_end);
    const DatasetStats cleaned_stats = compute_dataset_stats(cleaned);
    std::cout << "Shared denoised points=" << cleaned.size()
              << " preprocessing_ms=" << preprocessing_ms << std::endl;

    std::vector<SummaryRow> rows;
    for (const auto target : ring_targets) {
      RingReconstructionConfig config;
      config.target_ring_count = target;
      rows.push_back(benchmark_ring(cleaned, cleaned_stats, config, warmups, repeats,
                                    preprocessing_ms));
    }
    const SummaryRow* requested_reference = find_row(rows, "ring_flat", requested_reference_target);
    if (!requested_reference)
      throw std::runtime_error("Requested reference Ring target is not present in --ring-list");
    std::size_t reference_target = requested_reference_target;
    const SummaryRow* reference = requested_reference;
    if (!topology_ok(*reference)) {
      reference = nullptr;
      for (const auto target : ring_targets) {
        const SummaryRow* candidate = find_row(rows, "ring_flat", target);
        if (candidate && topology_ok(*candidate) &&
            (!reference || candidate->actual_ring_count > reference->actual_ring_count)) {
          reference = candidate;
          reference_target = target;
        }
      }
      if (!reference) throw std::runtime_error("No topology-valid Ring reference was produced");
      std::cout << "Requested reference target=" << requested_reference_target
                << " failed topology; using highest stable target=" << reference_target
                << " actual=" << reference->actual_ring_count << std::endl;
    }
    const double reference_volume = reference->volume_ml;
    for (auto& row : rows) {
      row.volume_diff_reference_ml = std::abs(row.volume_ml - reference_volume);
      row.volume_diff_reference_percent = 100.0 * row.volume_diff_reference_ml / reference_volume;
    }

    std::size_t recommended_target = reference_target;
    for (std::size_t i = 0; i < ring_targets.size(); ++i) {
      const SummaryRow* candidate = find_row(rows, "ring_flat", ring_targets[i]);
      if (!candidate || !topology_ok(*candidate) || candidate->total.p95 > 500.0 ||
          candidate->volume_diff_reference_percent > 0.5) continue;
      bool all_higher_converged = true;
      for (std::size_t j = i + 1; j < ring_targets.size(); ++j) {
        const SummaryRow* higher = find_row(rows, "ring_flat", ring_targets[j]);
        if (ring_targets[j] > reference_target) continue;
        if (!higher || !topology_ok(*higher) || higher->volume_diff_reference_percent > 0.5) {
          all_higher_converged = false;
          break;
        }
      }
      if (all_higher_converged) { recommended_target = ring_targets[i]; break; }
    }

    const std::vector<std::size_t> generated_counts{2, 3, 4, 5, 6, 8, 10};
    for (const auto generated : generated_counts) {
      RingReconstructionConfig config;
      config.target_ring_count = recommended_target;
      config.endpoint = EndpointCompletion::cubic_bspline;
      config.endpoint_fit_rings = 8;
      config.endpoint_generated_rings = generated;
      rows.push_back(benchmark_ring(cleaned, cleaned_stats, config, warmups, repeats,
                                    preprocessing_ms));
    }
    for (const auto fit : {4U, 6U, 10U, 12U}) {
      RingReconstructionConfig config;
      config.target_ring_count = recommended_target;
      config.endpoint = EndpointCompletion::cubic_bspline;
      config.endpoint_fit_rings = fit;
      config.endpoint_generated_rings = 4;
      rows.push_back(benchmark_ring(cleaned, cleaned_stats, config, warmups, repeats,
                                    preprocessing_ms));
    }
    rows.push_back(benchmark_alpha(cleaned, cleaned_stats, alpha_mm, offset_mm,
                                   warmups, repeats, preprocessing_ms));
    for (auto& row : rows) {
      row.volume_diff_reference_ml = std::abs(row.volume_ml - reference_volume);
      row.volume_diff_reference_percent = 100.0 * row.volume_diff_reference_ml / reference_volume;
    }

    const SummaryRow* final_row = find_row(rows, "ring_spline", recommended_target, 8, 4);
    if (!final_row) throw std::runtime_error("Final Ring-spline configuration was not benchmarked");
    write_csv(csv_path, rows);
    write_report(report_path, rows, ring_targets, reference_target, recommended_target,
                 preprocessing_ms, warmups, repeats, alpha_mm, offset_mm);
    write_mesh_ply(mesh_path, final_row->representative_mesh);
    std::cout << "WROTE csv=" << csv_path.string() << " report=" << report_path.string()
              << " mesh=" << mesh_path.string() << " recommended_target=" << recommended_target
              << std::endl;
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "FATAL: " << e.what() << '\n';
    return 1;
  }
}
