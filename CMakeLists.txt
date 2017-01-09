cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_INSTALL_PREFIX )

project(Bach)

# Setup ROOT as a dependency
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Physics MathCore Geom GenVector)

include(${ROOT_USE_FILE})

add_subdirectory(Millepede)
add_subdirectory(TbKernel)
add_subdirectory(TbAlgorithms)
add_subdirectory(Bach)
