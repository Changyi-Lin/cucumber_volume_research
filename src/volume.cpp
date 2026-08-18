#include "volume.hpp"

#include <CGAL/boost/graph/iterator.h>

#include <cmath>

std::pair<double, double> mesh_volume(const Mesh& mesh) {
  long double six_v = 0.0L;
  for (const auto f : mesh.faces()) {
    const auto h0 = mesh.halfedge(f);
    const auto h1 = mesh.next(h0);
    const auto h2 = mesh.next(h1);
    const Point& a = mesh.point(mesh.target(h0));
    const Point& b = mesh.point(mesh.target(h1));
    const Point& c = mesh.point(mesh.target(h2));
    const long double ax = a.x(), ay = a.y(), az = a.z();
    const long double bx = b.x(), by = b.y(), bz = b.z();
    const long double cx = c.x(), cy = c.y(), cz = c.z();
    six_v += ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) +
             az * (bx * cy - by * cx);
  }
  const double signed_v = static_cast<double>(six_v / 6.0L);
  return {signed_v, std::abs(signed_v)};
}

