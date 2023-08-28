cmake_minimum_required(VERSION 3.1...3.23)
project(Circular_kernel_3_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

find_package(Qt6 QUIET COMPONENTS Widgets OpenGL)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)

  add_executable(Circular_kernel_3 Circular_kernel_3.cpp Viewer.cpp)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS Circular_kernel_3)

  target_link_libraries(Circular_kernel_3 PRIVATE CGAL::CGAL CGAL::CGAL_Qt6
                                                  Qt6::Widgets Qt6::OpenGL)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Circular_kernel_3)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Circular_kernel_3)

else()

  message("NOTICE: This demo requires CGAL and Qt6, and will not be compiled.")

endif()
