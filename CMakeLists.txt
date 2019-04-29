cmake_minimum_required(VERSION 2.8.12)

project(inv_dist LANGUAGES CXX)

add_executable(inv_dist src/main.cpp)

set_target_properties(inv_dist PROPERTIES PREFIX "" CACHE STRING "prefix" FORCE)
set_target_properties(inv_dist PROPERTIES SUFFIX ".exe" CACHE STRING "suffix" FORCE)

include_directories(src)
#include_directories(lib/eigen)

SET ( CMAKE_BUILD_TYPE "Release" CACHE STRING "build_type" FORCE)
SET ( CMAKE_CXX_FLAGS "-std=c++17 -fopenmp -lpthread" CACHE STRING "cxx_flags" FORCE)
SET ( CMAKE_EXE_LINKER_FLAGS "-fopenmp -lpthread" CACHE STRING "linker_flags" FORCE)

set_target_properties(inv_dist
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ./bin
    LIBRARY_OUTPUT_DIRECTORY ./bin
    RUNTIME_OUTPUT_DIRECTORY ./bin
)
