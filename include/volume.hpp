#pragma once

#include "types.hpp"

#include <utility>

// Returns signed and absolute volume in mm^3.
std::pair<double, double> mesh_volume(const Mesh& mesh);

