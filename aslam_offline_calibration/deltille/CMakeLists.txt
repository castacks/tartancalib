## Copyright (C) 2017-present, Facebook, Inc.
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
##
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

project(Deltille VERSION 1.0.0 LANGUAGES C CXX)

# set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake")

# include(GNUInstallDirs)

# option(BUILD_APPS "Build the executable to detect and write corners to disk. Requires Boost" OFF)
# option(WITH_MATLAB "Build Matlab buindings" OFF)

# Regarding build type
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
else()
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "CMAKE_BUILD_TYPE is set to Debug. This will be produce inefficient code")
  endif()
endif()

set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall)

if(BUILD_APPS)
 add_subdirectory(apps)
endif()

# if(WITH_MATLAB)
#   include(ExternalProject)
#   add_subdirectory(matlab)
# endif()
