message("XSecVary Library")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -funroll-loops -march=native")

# Set the default verbosity of the compilation...
set(CMAKE_VERBOSE_MAKEFILE TRUE)

###################################################

# Configure the dependencies
find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

# Define the source and include files that should be used for the io
# part of ERepSim.

set(source_xsecvary
  XsecVary.cc XsecVary2D.cc SRGlobal.cc
  )

set(includes_xsecvary
  XsecVary.h XsecVary2D.h SRGlobal.h
  )

add_library(XsecVary SHARED ${source_xsecvary})

target_include_directories(XsecVary PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:include/XsecVary>")

install(TARGETS XsecVary
  EXPORT XsecVaryTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include/XsecVary)

# Install the header files.
install(FILES ${includes_XsecVary} DESTINATION include/XsecVary)


# Make sure the current directories are available for the later
# compilation.  This is required to make ROOT_GENERATE_DICTIONARY
# work.
# message(CMAKE_LIBRARY_OUTPUT_DIRECTORY="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

ROOT_GENERATE_DICTIONARY(SRGlobal_dict SRGlobal.h MODULE XsecVary LINKDEF LinkDef.h)
# install(FILES libXsecVary_rdict.pcm libXsecVary_rdict.rootmap DESTINATION "${CMAKE_BINARY_DIR}/lib" COMPONENT libraries)
install(TARGETS SRGlobal_dict DESTINATION lib COMPONENT libraries)
