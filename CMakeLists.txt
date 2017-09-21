cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_INSTALL_PREFIX )

project(Bach)

include ( ${DD4hep_DIR}/cmake/DD4hep.cmake )

dd4hep_configure_output()

# Setup ROOT as a dependency
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Physics MathCore Geom GenVector)
find_package(DD4hep REQUIRED COMPONENTS DDCore DDCond DDAlign)

include(${ROOT_USE_FILE})

add_subdirectory(Millepede)
add_subdirectory(TbKernel)
add_subdirectory(TbAlgorithms)
add_subdirectory(Bach)
