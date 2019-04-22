cmake_minimum_required(VERSION 3.1.0)

project(HilbertCurve CXX)

add_compile_options(-Wall -Wextra -Werror)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  add_compile_options(-Ofast -march=native -mtune=native)
endif()

add_executable(example example.cpp)
add_executable(hilbert_test hilbert_test.cpp)

set_target_properties(example PROPERTIES CXX_STANDARD 17)
set_target_properties(hilbert_test PROPERTIES CXX_STANDARD 17)
