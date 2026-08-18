#pragma once

#include "types.hpp"

#include <vector>

Mesh run_alpha_wrap(const std::vector<Point>& points, double alpha_mm,
                    double offset_mm);

