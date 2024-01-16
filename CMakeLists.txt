cmake_minimum_required(VERSION 3.0)

# Set the project name and version.  The project name will be used to
# import this into other CMAKE files.
project(XsecResponse VERSION 1.0.0)

message("XsecResponse -- ${VERSION}")

# Define the options that can be set in the cache, or on the cmake
# command line.
set(CMAKE_BUILD_TYPE Debug)

# Changes default install path to be a subdirectory of the build dir.
# Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")

# Make sure that ROOT is available.  ROOT is absolutely required.
find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO EG)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

add_subdirectory(src_bin)
add_subdirectory(src_lib)