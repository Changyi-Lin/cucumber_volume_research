cmake_minimum_required(VERSION 3.20)
project(cucumber_alpha_wrap LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(TBB QUIET)
find_package(OpenMP QUIET)

add_executable(cucumber_alpha_wrap
  src/main.cpp
  src/ply_reader.cpp
  src/voxel_downsample.cpp
  src/alpha_wrap.cpp
  src/mesh_validation.cpp
  src/volume.cpp
  src/denoise.cpp
)

target_include_directories(cucumber_alpha_wrap PRIVATE include)
target_link_libraries(cucumber_alpha_wrap PRIVATE CGAL::CGAL Eigen3::Eigen)
if(TBB_FOUND)
  target_link_libraries(cucumber_alpha_wrap PRIVATE TBB::tbb)
  target_compile_definitions(cucumber_alpha_wrap PRIVATE CUCUMBER_HAS_TBB=1)
endif()
if(OpenMP_CXX_FOUND)
  target_link_libraries(cucumber_alpha_wrap PRIVATE OpenMP::OpenMP_CXX)
  target_compile_definitions(cucumber_alpha_wrap PRIVATE CUCUMBER_HAS_OPENMP=1)
endif()

if(MSVC)
  target_compile_options(cucumber_alpha_wrap PRIVATE /O2 /DNDEBUG)
else()
  target_compile_options(cucumber_alpha_wrap PRIVATE -O3 -march=native -DNDEBUG -Wall -Wextra -Wpedantic)
endif()
