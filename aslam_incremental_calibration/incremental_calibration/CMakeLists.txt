cmake_minimum_required(VERSION 3.0.2)
project(incremental_calibration)

set(CMAKE_CXX_STANDARD 14)

find_package(catkin REQUIRED COMPONENTS sm_eigen aslam_backend)
find_package(Boost REQUIRED COMPONENTS system thread)
find_package(Eigen3 REQUIRED)

catkin_package(
  LIBRARIES ${PROJECT_NAME}
  INCLUDE_DIRS include
  CATKIN_DEPENDS sm_eigen
)

include_directories(include ${catkin_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

add_library(${PROJECT_NAME}
  src/base/Serializable.cpp
  src/base/Timestamp.cpp
  src/base/Condition.cpp
  src/base/Mutex.cpp
  src/base/Thread.cpp
  src/base/Threads.cpp
  src/base/Timer.cpp
  src/exceptions/Exception.cpp
  src/exceptions/InvalidOperationException.cpp
  src/exceptions/NullPointerException.cpp
  src/exceptions/SystemException.cpp
  src/statistics/NormalDistribution1v.cpp
  src/statistics/ChiSquareDistribution.cpp
  src/statistics/EstimatorMLNormal1v.cpp
  src/functions/IncompleteGammaPFunction.cpp
  src/functions/IncompleteGammaQFunction.cpp
  src/functions/LogFactorialFunction.cpp
  src/functions/LogGammaFunction.cpp
  src/core/IncrementalEstimator.cpp
  src/core/OptimizationProblem.cpp
  src/core/IncrementalOptimizationProblem.cpp
  src/core/LinearSolver.cpp
  src/core/LinearSolverOptions.cpp
  src/algorithms/marginalize.cpp
  src/algorithms/linalg.cpp
)
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES} ${Boost_LIBRARIES} ${TBB_LIBRARIES})

# if(CATKIN_ENABLE_TESTING)
#   # Avoid clash with tr1::tuple:
#   # https://code.google.com/p/googletest/source/browse/trunk/README?r=589#257
#   add_definitions(-DGTEST_USE_OWN_TR1_TUPLE=0)

#   catkin_add_gtest(${PROJECT_NAME}_test
#     test/test_main.cpp
#     test/VectorDesignVariableTest.cpp
#     test/AlgorithmsTest.cpp
#     test/OptimizationProblemTest.cpp
#     test/IncrementalOptimizationProblemTest.cpp
#     test/LinearSolverTest.cpp
#     )
#   target_link_libraries(${PROJECT_NAME}_test ${PROJECT_NAME})
# endif()

##################
## Installation ##
##################

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

## Mark cpp header files for installation
install(DIRECTORY include/
  DESTINATION ${CATKIN_GLOBAL_INCLUDE_DESTINATION}
)
