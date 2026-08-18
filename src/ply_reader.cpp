#include "ply_reader.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

enum class PlyFormat { ascii, binary_little, binary_big };

struct Property {
  std::string type;
  std::string name;
  bool is_list = false;
};

struct Header {
  PlyFormat format = PlyFormat::ascii;
  std::size_t vertices = 0;
  std::vector<Property> vertex_properties;
};

std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return s;
}

Header parse_header(std::istream& in) {
  Header header;
  std::string line;
  if (!std::getline(in, line) || lower(line.substr(0, 3)) != "ply")
    throw std::runtime_error("Not a PLY file");
  std::string element;
  bool got_format = false;
  bool got_end = false;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    std::istringstream ss(line);
    std::string keyword;
    ss >> keyword;
    keyword = lower(keyword);
    if (keyword == "format") {
      std::string value;
      ss >> value;
      value = lower(value);
      if (value == "ascii") header.format = PlyFormat::ascii;
      else if (value == "binary_little_endian") header.format = PlyFormat::binary_little;
      else if (value == "binary_big_endian") header.format = PlyFormat::binary_big;
      else throw std::runtime_error("Unsupported PLY format: " + value);
      got_format = true;
    } else if (keyword == "element") {
      std::size_t count = 0;
      ss >> element >> count;
      element = lower(element);
      if (element == "vertex") header.vertices = count;
    } else if (keyword == "property" && element == "vertex") {
      Property p;
      ss >> p.type;
      p.type = lower(p.type);
      if (p.type == "list") {
        p.is_list = true;
        std::string count_type;
        ss >> count_type >> p.type >> p.name;
      } else {
        ss >> p.name;
      }
      p.name = lower(p.name);
      header.vertex_properties.push_back(std::move(p));
    } else if (keyword == "end_header") {
      got_end = true;
      break;
    }
  }
  if (!got_format || !got_end || header.vertices == 0)
    throw std::runtime_error("Incomplete PLY header");
  if (header.vertex_properties.empty())
    throw std::runtime_error("PLY vertex element has no properties");
  for (const auto& p : header.vertex_properties)
    if (p.is_list) throw std::runtime_error("List property in PLY vertex element is unsupported");
  return header;
}

template <class T>
T byte_swap(T value) {
  std::array<std::byte, sizeof(T)> bytes{};
  std::memcpy(bytes.data(), &value, sizeof(T));
  std::reverse(bytes.begin(), bytes.end());
  std::memcpy(&value, bytes.data(), sizeof(T));
  return value;
}

template <class T>
double read_binary_as_double(std::istream& in, const bool swap) {
  T value{};
  in.read(reinterpret_cast<char*>(&value), sizeof(value));
  if (!in) throw std::runtime_error("Unexpected EOF in binary PLY vertices");
  if (swap && sizeof(T) > 1) value = byte_swap(value);
  return static_cast<double>(value);
}

double read_scalar(std::istream& in, std::string type, const bool swap) {
  type = lower(std::move(type));
  if (type == "char" || type == "int8") return read_binary_as_double<std::int8_t>(in, false);
  if (type == "uchar" || type == "uint8") return read_binary_as_double<std::uint8_t>(in, false);
  if (type == "short" || type == "int16") return read_binary_as_double<std::int16_t>(in, swap);
  if (type == "ushort" || type == "uint16") return read_binary_as_double<std::uint16_t>(in, swap);
  if (type == "int" || type == "int32") return read_binary_as_double<std::int32_t>(in, swap);
  if (type == "uint" || type == "uint32") return read_binary_as_double<std::uint32_t>(in, swap);
  if (type == "float" || type == "float32") return read_binary_as_double<float>(in, swap);
  if (type == "double" || type == "float64") return read_binary_as_double<double>(in, swap);
  throw std::runtime_error("Unsupported PLY scalar type: " + type);
}

}  // namespace

std::vector<Point> read_ply_xyz(const std::filesystem::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("Cannot open PLY: " + path.string());
  const Header header = parse_header(in);
  int xi = -1, yi = -1, zi = -1;
  for (std::size_t i = 0; i < header.vertex_properties.size(); ++i) {
    if (header.vertex_properties[i].name == "x") xi = static_cast<int>(i);
    if (header.vertex_properties[i].name == "y") yi = static_cast<int>(i);
    if (header.vertex_properties[i].name == "z") zi = static_cast<int>(i);
  }
  if (xi < 0 || yi < 0 || zi < 0) throw std::runtime_error("PLY lacks x/y/z vertex properties");

  std::vector<Point> points;
  points.reserve(header.vertices);
  if (header.format == PlyFormat::ascii) {
    for (std::size_t row = 0; row < header.vertices; ++row) {
      double x = 0, y = 0, z = 0;
      for (std::size_t col = 0; col < header.vertex_properties.size(); ++col) {
        double value = 0;
        if (!(in >> value)) throw std::runtime_error("Unexpected EOF in ASCII PLY vertices");
        if (static_cast<int>(col) == xi) x = value;
        if (static_cast<int>(col) == yi) y = value;
        if (static_cast<int>(col) == zi) z = value;
      }
      points.emplace_back(x, y, z);
    }
  } else {
    const bool file_little = header.format == PlyFormat::binary_little;
    const std::uint16_t endian_probe = 1;
    const bool host_little = *reinterpret_cast<const std::uint8_t*>(&endian_probe) == 1;
    const bool swap = file_little != host_little;
    for (std::size_t row = 0; row < header.vertices; ++row) {
      double x = 0, y = 0, z = 0;
      for (std::size_t col = 0; col < header.vertex_properties.size(); ++col) {
        const double value = read_scalar(in, header.vertex_properties[col].type, swap);
        if (static_cast<int>(col) == xi) x = value;
        if (static_cast<int>(col) == yi) y = value;
        if (static_cast<int>(col) == zi) z = value;
      }
      points.emplace_back(x, y, z);
    }
  }
  return points;
}

DatasetStats compute_dataset_stats(const std::vector<Point>& points) {
  if (points.empty()) throw std::runtime_error("Cannot compute statistics of empty point set");
  DatasetStats s;
  s.point_count = points.size();
  s.min = {{points[0].x(), points[0].y(), points[0].z()}};
  s.max = s.min;
  Eigen::Vector3d sum = Eigen::Vector3d::Zero();
  for (const auto& p : points) {
    const Eigen::Vector3d v(p.x(), p.y(), p.z());
    sum += v;
    for (int k = 0; k < 3; ++k) {
      s.min[k] = std::min(s.min[k], v[k]);
      s.max[k] = std::max(s.max[k], v[k]);
    }
  }
  const Eigen::Vector3d mean = sum / static_cast<double>(points.size());
  for (int k = 0; k < 3; ++k) s.centroid[k] = mean[k];

  Eigen::Matrix3d covariance = Eigen::Matrix3d::Zero();
  for (const auto& p : points) {
    const Eigen::Vector3d d(p.x() - mean.x(), p.y() - mean.y(), p.z() - mean.z());
    covariance.noalias() += d * d.transpose();
  }
  covariance /= static_cast<double>(points.size());
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(covariance);
  if (solver.info() != Eigen::Success) throw std::runtime_error("PCA eigensolver failed");
  const Eigen::Matrix3d axes = solver.eigenvectors();
  for (int col = 0; col < 3; ++col) {
    s.eigenvalues[col] = solver.eigenvalues()[col];
    for (int row = 0; row < 3; ++row) s.pca_axes[col][row] = axes(row, col);
  }
  Eigen::Vector3d pmin = Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());
  Eigen::Vector3d pmax = Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity());
  for (const auto& p : points) {
    const Eigen::Vector3d d(p.x() - mean.x(), p.y() - mean.y(), p.z() - mean.z());
    const Eigen::Vector3d q = axes.transpose() * d;
    pmin = pmin.cwiseMin(q);
    pmax = pmax.cwiseMax(q);
  }
  for (int k = 0; k < 3; ++k) s.pca_extents[k] = pmax[k] - pmin[k];
  const double dx = s.max[0] - s.min[0], dy = s.max[1] - s.min[1], dz = s.max[2] - s.min[2];
  s.aabb_diagonal = std::sqrt(dx * dx + dy * dy + dz * dz);
  return s;
}
