# hilbert-curve: Optimized N-dimensional Hilbert curve generation in C++
# Copyright (C) 2019 <tomKPZ@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

cmake_minimum_required(VERSION 3.1.0)

project(HilbertCurve CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall -Wextra -Werror)
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  add_compile_options(-Ofast -march=native -mtune=native)
endif()

add_executable(benchmark benchmark.cpp)
add_executable(example example.cpp)
add_executable(hilbert_test hilbert_test.cpp)
add_executable(temp temp.cpp)
