#include "alpha_wrap.hpp"

#include <CGAL/alpha_wrap_3.h>

#include <stdexcept>

Mesh run_alpha_wrap(const std::vector<Point>& points, const double alpha_mm,
                    const double offset_mm) {
  if (points.empty()) throw std::runtime_error("Alpha Wrap input is empty");
  if (!(alpha_mm > 0) || !(offset_mm > 0))
    throw std::runtime_error("alpha and offset must be positive");
  Mesh output;
  CGAL::alpha_wrap_3(points, alpha_mm, offset_mm, output);
  return output;
}

