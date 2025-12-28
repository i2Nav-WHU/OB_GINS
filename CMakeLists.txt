cmake_minimum_required(VERSION 3.10)

set(CMAKE_TOOLCHAIN_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake"
    CACHE STRING "Vcpkg toolchain file")

project(ob_gins)

set(CMAKE_CXX_STANDARD 17)
if(MSVC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /utf-8")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /utf-8")
  add_compile_options(/permissive-)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

add_compile_definitions(M_PI=3.14159265358979323846)

file(
  GLOB_RECURSE
  SOURCE
  src/ob_gins.cc
  src/fileio/fileloader.cc
  src/fileio/filesaver.cc
  src/factors/pose_manifold.cc
  src/preintegration/preintegration_base.cc
  src/preintegration/preintegration_normal.cc
  src/preintegration/preintegration_earth.cc
  src/preintegration/preintegration_odo.cc
  src/preintegration/preintegration_earth_odo.cc)

include_directories(${PROJECT_SOURCE_DIR})
add_executable(${PROJECT_NAME} ${SOURCE})

# Eigen3
find_package(Eigen3 REQUIRED)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

# yaml-cpp
find_package(yaml-cpp REQUIRED)
target_link_libraries(${PROJECT_NAME} yaml-cpp::yaml-cpp)

# abseil
set(ABSL_PROPAGATE_CXX_STD true)
add_subdirectory(src/thirdparty/abseil-cpp abseil-cpp)
target_link_libraries(${PROJECT_NAME} absl::strings absl::str_format absl::time)

# Ceres
find_package(Ceres REQUIRED)
target_link_libraries(${PROJECT_NAME} Ceres::ceres)
