cmake_minimum_required(VERSION 3.0.2)
project(tartan_calib)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 REQUIRED NO_MODULE)
if(Eigen3_FOUND)
    message(STATUS "EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}")
endif()

## libvis





find_package(GLEW REQUIRED)
set( CMAKE_CXX_COMPILE_OPTIONS_PIE "-fPIC" )


find_package(Boost REQUIRED COMPONENTS serialization) 
find_package(OpenMP REQUIRED)
find_package(catkin REQUIRED COMPONENTS
  sm_common
  sm_eigen
  sm_boost
  sm_kinematics
  aslam_time
  sm_opencv
  sm_logging
  opencv2_catkin
  aslam_cameras
  sm_property_tree
  aslam_time
  ethz_apriltag2
  aslam_cameras_april
)

catkin_package(
  INCLUDE_DIRS include
)

add_definitions(-fPIC -Wextra -Winit-self -Woverloaded-virtual -Wnon-virtual-dtor -Wsign-promo -Wno-long-long -Wl,--no-allow-shlib-undefined)

set(OTHER_LIBS)
if(NOT APPLE)
  set(OTHER_LIBS rt)
endif()



add_library(tartan_calib INTERFACE)



