#include "mesh_validation.hpp"

#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/boost/graph/helpers.h>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <queue>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;

namespace {

std::vector<std::size_t> component_face_counts(const Mesh& mesh) {
  if (mesh.number_of_faces() == 0) return {};
  std::size_t max_index = 0;
  for (const auto face : mesh.faces()) max_index = std::max(max_index, static_cast<std::size_t>(face.idx()));
  std::vector<unsigned char> seen(max_index + 1, 0);
  std::vector<std::size_t> counts;
  std::queue<Mesh::Face_index> queue;
  for (const auto start : mesh.faces()) {
    if (seen[start.idx()]) continue;
    counts.push_back(0);
    seen[start.idx()] = 1;
    queue.push(start);
    while (!queue.empty()) {
      const auto face = queue.front();
      queue.pop();
      ++counts.back();
      for (const auto h : CGAL::halfedges_around_face(mesh.halfedge(face), mesh)) {
        const auto opposite_face = mesh.face(mesh.opposite(h));
        if (opposite_face != Mesh::null_face() && !seen[opposite_face.idx()]) {
          seen[opposite_face.idx()] = 1;
          queue.push(opposite_face);
        }
      }
    }
  }
  std::sort(counts.begin(), counts.end(), std::greater<std::size_t>());
  return counts;
}

double triangle_area(const Point& a, const Point& b, const Point& c) {
  const Eigen::Vector3d ab(b.x() - a.x(), b.y() - a.y(), b.z() - a.z());
  const Eigen::Vector3d ac(c.x() - a.x(), c.y() - a.y(), c.z() - a.z());
  return 0.5 * ab.cross(ac).norm();
}

bool line_triangle_intersection(const Eigen::Vector3d& origin,
                                const Eigen::Vector3d& direction,
                                const Point& pa, const Point& pb, const Point& pc,
                                double& t) {
  const Eigen::Vector3d a(pa.x(), pa.y(), pa.z());
  const Eigen::Vector3d b(pb.x(), pb.y(), pb.z());
  const Eigen::Vector3d c(pc.x(), pc.y(), pc.z());
  const Eigen::Vector3d e1 = b - a;
  const Eigen::Vector3d e2 = c - a;
  const Eigen::Vector3d p = direction.cross(e2);
  const double det = e1.dot(p);
  if (std::abs(det) < 1e-12) return false;
  const double inv_det = 1.0 / det;
  const Eigen::Vector3d s = origin - a;
  const double u = s.dot(p) * inv_det;
  if (u < -1e-10 || u > 1.0 + 1e-10) return false;
  const Eigen::Vector3d q = s.cross(e1);
  const double v = direction.dot(q) * inv_det;
  if (v < -1e-10 || u + v > 1.0 + 1e-10) return false;
  t = e2.dot(q) * inv_det;
  return std::isfinite(t);
}

std::size_t unique_line_intersections(const Mesh& mesh,
                                      const Eigen::Vector3d& origin,
                                      const Eigen::Vector3d& direction,
                                      const double cluster_tolerance) {
  std::vector<double> hits;
  hits.reserve(16);
  for (const auto face : mesh.faces()) {
    const auto h0 = mesh.halfedge(face);
    const auto h1 = mesh.next(h0);
    const auto h2 = mesh.next(h1);
    double t = 0;
    if (line_triangle_intersection(origin, direction,
                                   mesh.point(mesh.target(h0)),
                                   mesh.point(mesh.target(h1)),
                                   mesh.point(mesh.target(h2)), t))
      hits.push_back(t);
  }
  if (hits.empty()) return 0;
  std::sort(hits.begin(), hits.end());
  std::size_t unique = 1;
  double previous = hits.front();
  for (std::size_t i = 1; i < hits.size(); ++i) {
    if (hits[i] - previous > cluster_tolerance) ++unique;
    previous = hits[i];
  }
  return unique;
}

void detect_two_sided(const Mesh& mesh, const std::vector<Point>& input,
                      const DatasetStats& stats, const double offset_mm,
                      ValidationResult& out) {
  if (mesh.number_of_faces() == 0 || input.empty()) return;
  Eigen::Vector3d mean(stats.centroid[0], stats.centroid[1], stats.centroid[2]);
  Eigen::Vector3d transverse0(stats.pca_axes[0][0], stats.pca_axes[0][1], stats.pca_axes[0][2]);
  Eigen::Vector3d transverse1(stats.pca_axes[1][0], stats.pca_axes[1][1], stats.pca_axes[1][2]);
  Eigen::Vector3d axis(stats.pca_axes[2][0], stats.pca_axes[2][1], stats.pca_axes[2][2]);
  transverse0.normalize(); transverse1.normalize(); axis.normalize();

  double axial_min = std::numeric_limits<double>::infinity();
  double axial_max = -std::numeric_limits<double>::infinity();
  for (const auto& p : input) {
    const Eigen::Vector3d v(p.x(), p.y(), p.z());
    const double s = axis.dot(v - mean);
    axial_min = std::min(axial_min, s);
    axial_max = std::max(axial_max, s);
  }
  const double span = axial_max - axial_min;
  const double slab_half_width = std::max(0.75, span * 0.0125);
  const double cluster_tolerance = std::max(1e-5, offset_mm * 0.002);
  constexpr double fractions[] = {0.20, 0.35, 0.50, 0.65, 0.80};
  std::vector<std::size_t> counts;
  for (const double fraction : fractions) {
    const double target = axial_min + fraction * span;
    Eigen::Vector3d local_sum = Eigen::Vector3d::Zero();
    std::size_t local_count = 0;
    for (const auto& p : input) {
      const Eigen::Vector3d v(p.x(), p.y(), p.z());
      if (std::abs(axis.dot(v - mean) - target) <= slab_half_width) {
        local_sum += v;
        ++local_count;
      }
    }
    if (local_count < 3) continue;
    const Eigen::Vector3d local_center = local_sum / static_cast<double>(local_count);
    counts.push_back(unique_line_intersections(mesh, local_center, transverse0, cluster_tolerance));
    counts.push_back(unique_line_intersections(mesh, local_center, transverse1, cluster_tolerance));
  }
  if (counts.empty()) return;
  out.tested_lines = counts.size();
  out.lines_with_four_or_more = static_cast<std::size_t>(
      std::count_if(counts.begin(), counts.end(), [](std::size_t n) { return n >= 4; }));
  std::sort(counts.begin(), counts.end());
  const std::size_t mid = counts.size() / 2;
  out.median_line_intersections = counts.size() % 2
      ? static_cast<double>(counts[mid])
      : 0.5 * static_cast<double>(counts[mid - 1] + counts[mid]);
  const double four_ratio = static_cast<double>(out.lines_with_four_or_more) / out.tested_lines;
  out.two_sided_wrap = out.tested_lines >= 6 &&
                       out.median_line_intersections >= 4.0 && four_ratio >= 0.4;
}

}  // namespace

ComponentCleanupResult cleanup_tiny_components(Mesh& mesh,
                                               const double minimum_largest_face_ratio) {
  ComponentCleanupResult result;
  const auto counts = component_face_counts(mesh);
  result.raw_connected_components = counts.size();
  if (counts.empty()) {
    result.largest_component_face_ratio = 0;
    return result;
  }
  result.largest_component_face_ratio =
      static_cast<double>(counts.front()) / static_cast<double>(mesh.number_of_faces());
  if (counts.size() > 1 && result.largest_component_face_ratio >= minimum_largest_face_ratio) {
    result.removed_components = PMP::keep_largest_connected_components(mesh, 1);
    mesh.collect_garbage();
  }
  return result;
}

ValidationResult validate_mesh(const Mesh& mesh, const std::vector<Point>& input,
                               const DatasetStats& input_stats,
                               const double offset_mm) {
  ValidationResult result;
  result.vertices = mesh.number_of_vertices();
  result.faces = mesh.number_of_faces();
  result.edges = mesh.number_of_edges();
  result.is_closed = CGAL::is_closed(mesh);
  result.is_triangle_mesh = CGAL::is_triangle_mesh(mesh);
  result.connected_components = component_face_counts(mesh).size();
  result.euler_characteristic = static_cast<long long>(result.vertices) -
                                static_cast<long long>(result.edges) +
                                static_cast<long long>(result.faces);
  if (result.is_closed && result.connected_components > 0) {
    result.genus = (2.0 * static_cast<double>(result.connected_components) -
                    static_cast<double>(result.euler_characteristic)) / 2.0;
  }

  std::vector<Mesh::Halfedge_index> non_manifold_vertices;
  if (CGAL::is_valid_polygon_mesh(mesh))
    PMP::non_manifold_vertices(mesh, std::back_inserter(non_manifold_vertices));
  result.manifold = CGAL::is_valid_polygon_mesh(mesh) && non_manifold_vertices.empty();

  if (result.is_triangle_mesh && result.faces > 0) {
    result.self_intersection = PMP::does_self_intersect(mesh);
    if (result.is_closed) result.outward_oriented = PMP::is_outward_oriented(mesh);
  }

  double area = 0;
  for (const auto face : mesh.faces()) {
    const auto h0 = mesh.halfedge(face);
    const auto h1 = mesh.next(h0);
    const auto h2 = mesh.next(h1);
    area += triangle_area(mesh.point(mesh.target(h0)),
                          mesh.point(mesh.target(h1)),
                          mesh.point(mesh.target(h2)));
  }
  result.surface_area_mm2 = area;
  detect_two_sided(mesh, input, input_stats, offset_mm, result);
  return result;
}
