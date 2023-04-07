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

project(deltille VERSION 1.0.0 LANGUAGES C CXX)

# set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake")

# include(GNUInstallDirs)

# option(BUILD_APPS "Build the executable to detect and write corners to disk. Requires Boost" OFF)
# option(WITH_MATLAB "Build Matlab buindings" OFF)

find_package(catkin)
find_package(OpenCV CONFIG REQUIRED COMPONENTS core highgui imgproc imgcodecs) 

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
)

# Regarding build type
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
else()
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "CMAKE_BUILD_TYPE is set to Debug. This will be produce inefficient code")
  endif()
endif()

## [START COPY] ##

add_library(opencv INTERFACE IMPORTED)
set_target_properties(opencv PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${OpenCV_INCLUDE_DIRS}"
  INTERFACE_LINK_LIBRARIES "${OpenCV_LIBS}")

add_library(${PROJECT_NAME}
  src/DetectorTools.cpp
  src/TagFamily.cpp
  src/TaggedBoardIndexer.cpp
  src/PolynomialFit.cpp)

target_include_directories(deltille PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(${PROJECT_NAME} PUBLIC ${catkin_LIBRARIES} opencv)
target_compile_options(${PROJECT_NAME} PUBLIC "-Wno-deprecated-copy") # TODO: Fix...
# target_compile_features(deltille PUBLIC cxx_auto_type)

install(TARGETS deltille EXPORT ${PROJECT_NAME}Config
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)

# install(DIRECTORY ${CMAKE_SOURCE_DIR}/include DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT "${PROJECT_NAME}Config" DESTINATION share/${PROJECT_NAME}/cmake)
export(TARGETS deltille FILE ${PROJECT_NAME}Config.cmake)

## [END COPY] ##

set(CMAKE_CXX_EXTENSIONS OFF)

add_compile_options(-Wall)

# if(BUILD_APPS)
#  add_subdirectory(apps)
# endif()

# if(WITH_MATLAB)
#   include(ExternalProject)
#   add_subdirectory(matlab)
# endif()
