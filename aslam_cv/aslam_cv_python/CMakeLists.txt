cmake_minimum_required(VERSION 3.0.2)
project(aslam_cv_python)

set(CMAKE_CXX_STANDARD 14)

# libvis


set(LIBVIS_FILES
  libvis/src/libvis/camera.h
  libvis/src/libvis/camera_frustum_opengl.h
  libvis/src/libvis/command_line_parser.cc
  libvis/src/libvis/command_line_parser.h
  libvis/src/libvis/dlt.h
  libvis/src/libvis/eigen.h
  libvis/src/libvis/geometry.h
  libvis/src/libvis/glew.cc
  libvis/src/libvis/glew.h
  libvis/src/libvis/image.cc
  libvis/src/libvis/image.h
  libvis/src/libvis/image_cache.h
  libvis/src/libvis/image_frame.h
  libvis/src/libvis/image_io.cc
  libvis/src/libvis/image_io.h
  libvis/src/libvis/image_io_libpng.cc
  libvis/src/libvis/image_io_libpng.h
  libvis/src/libvis/image_io_netpbm.cc
  libvis/src/libvis/image_io_netpbm.h
  libvis/src/libvis/lm_optimizer.h
  libvis/src/libvis/lm_optimizer_impl.h
  libvis/src/libvis/lm_optimizer_update_accumulator.h
  libvis/src/libvis/loss_functions.h
  libvis/src/libvis/mesh.h
  libvis/src/libvis/mesh_opengl.h
  libvis/src/libvis/opengl.cc
  libvis/src/libvis/opengl.h
  libvis/src/libvis/opengl_context_glx.cc
  libvis/src/libvis/opengl_context_glx.h
  libvis/src/libvis/point_cloud.h
  libvis/src/libvis/point_cloud_opengl.h
  libvis/src/libvis/rgbd_video.h
  libvis/src/libvis/rgbd_video_io_tum_dataset.h
  libvis/src/libvis/shader_program_opengl.cc
  libvis/src/libvis/shader_program_opengl.h
  libvis/src/libvis/statistics.h
  libvis/src/libvis/timing.cc
  libvis/src/libvis/timing.h
  libvis/src/libvis/util.h
  
  ${GENERATED_HEADERS}
  libvis/resources/resources.qrc
)

add_library(libvis SHARED
  ${LIBVIS_FILES}
)
set_property(TARGET libvis PROPERTY FOLDER "libvis")


set(DELTILLE_FILES 
deltille/deltille.cc  
deltille/DetectorParams.h 
deltille/PolynomialFit.h 
deltille/TagFamily.h
deltille/utils.h
deltille/PolynomialSaddleDetectorContext.h
deltille/DetectorTools.cpp
)

add_library(deltille SHARED IMPORTED
${DELTILLE_FILES}
)

# set_property(TARGET deltille PROPERTY FOLDER "deltille")

# set_property(TARGET libvis PROPERTY FOLDER "libvis")
include_directories(
  ./deltille/
  ./libvis/src
  ./libvis/third_party/loguru/
    include 
    ${catkin_INCLUDE_DIRS}
    ${GLEW_INCLUDE_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${Qt5Core_INCLUDE_DIRS}
)

# add_library(loguru SHARED
#   libvis/third_party/loguru/loguru.cpp
# )
# set_property(TARGET loguru PROPERTY FOLDER "libvis/third_party")
# include_directories(libvis/third_party/loguru/)

# set_target_properties(loguru PROPERTIES
#   POSITION_INDEPENDENT_CODE ON
# )
# target_compile_definitions(loguru PRIVATE
#   LOGURU_WITH_STREAMS
# )

find_package(Boost REQUIRED COMPONENTS serialization) 
find_package(catkin REQUIRED COMPONENTS
  ethz_apriltag2  
  tartan_calib
  sm_common
  sm_eigen
  sm_kinematics
  sm_python
  sm_opencv
  opencv2_catkin
  python_module
  aslam_cv_serialization
  aslam_imgproc
  aslam_cameras
  aslam_time
  aslam_cv_serialization
  aslam_cameras_april
)

catkin_package(
  INCLUDE_DIRS include
)

add_definitions(-fPIC -Wextra -Winit-self -Woverloaded-virtual -Wnon-virtual-dtor -Wsign-promo -Wno-long-long -Wl,--no-allow-shlib-undefined)

if(APPLE)
  add_definitions( -ftemplate-depth-1024)
endif()

include_directories(include ${catkin_INCLUDE_DIRS})

add_python_export_library(${PROJECT_NAME} python/aslam_cv
  src/module.cpp
  src/tartancalibworker_.cpp
  src/CameraGeometries.cpp
  src/CameraGeometry.cpp
  src/Time.cpp
#  src/CameraSystem.cpp
#  src/CameraSystemClasses.cpp
#  src/MatchingAlgorithms.cpp
#  src/MultiFrame.cpp
  src/Frontend.cpp
  src/Frame.cpp
#  src/Descriptors.cpp
  src/CameraProjections.cpp
  src/CameraShutters.cpp
#  src/ImageSynchronizer.cpp
#  src/FrameBuilder.cpp
#  src/ImageContainer.cpp
  src/GridCalibration.cpp
#  src/Landmark.cpp
  src/Undistorters.cpp
#  src/NCameras.cpp
#  src/MultiFrameFactory.cpp
  src/OmniUndistorter.cpp
  src/PinholeUndistorter.cpp
)
target_compile_options(${PROJECT_NAME} PUBLIC "-Wno-deprecated-copy") # TODO: Fix...
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

##################
## Installation ##
##################
install(DIRECTORY include/
  DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION}
)