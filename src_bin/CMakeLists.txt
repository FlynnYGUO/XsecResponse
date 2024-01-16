# Build all of the applications here.
project(app)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -funroll-loops -march=native")

# Set the default verbosity of the compilation...
set(CMAKE_VERBOSE_MAKEFILE TRUE)


find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

add_executable(make_xsec_response_1d make_xsec_response_1d.cc)
target_link_libraries(make_xsec_response_1d
  LINK_PUBLIC XsecVary ROOT::EG ROOT::RIO ${ROOT_LIBRARIES})

add_executable(make_xsec_response_2d make_xsec_response_2d.cc)
target_link_libraries(make_xsec_response_2d
LINK_PUBLIC XsecVary ROOT::EG ROOT::RIO ${ROOT_LIBRARIES})


install(TARGETS make_xsec_response_1d make_xsec_response_2d RUNTIME DESTINATION bin)

