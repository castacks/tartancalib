ExternalProject_Add(MexMat
  GIT_REPOSITORY "https://github.com/halismai/mexmat.git"
  PREFIX         "${CMAKE_CURRENT_BINARY_DIR}/mexmat"
  PATCH_COMMAND   ""
  INSTALL_COMMAND ""
  BUILD_COMMAND   ""
  CONFIGURE_COMMAND ""
  SOURCE_DIR     "${CMAKE_CURRENT_BINARY_DIR}/mexmat")

set(MexMat_INCLUDE_DIRS "${CMAKE_CURRENT_BINARY_DIR}")
add_library(mexmat INTERFACE)

target_include_directories(mexmat INTERFACE
  $<BUILD_INTERFACE:${MexMat_INCLUDE_DIRS}>
  $<INSTALL_INTERFACE:include>)

target_compile_definitions(mexmat INTERFACE 
  -DMEXMAT_WITH_OPENCV)



