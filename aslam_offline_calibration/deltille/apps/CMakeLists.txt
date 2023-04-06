find_package(Boost REQUIRED COMPONENTS filesystem program_options system)

add_library(boost INTERFACE IMPORTED)
set_target_properties(boost PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
  INTERFACE_LINK_LIBRARIES "${Boost_LIBRARIES}")

add_executable(deltille_detector DeltilleDetector.cpp)
target_link_libraries(deltille_detector deltille boost)

install(TARGETS deltille_detector DESTINATION ${CMAKE_INSTALL_BINDIR})
