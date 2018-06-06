project (Circular_kernel_3_Demo)

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

find_package(CGAL COMPONENTS Qt5)
include(${CGAL_USE_FILE})

find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND)


  include_directories (BEFORE ../../include ./ )

  add_executable (Circular_kernel_3 Circular_kernel_3.cpp Viewer.cpp ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES})

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Circular_kernel_3 )

  target_link_libraries( Circular_kernel_3 PRIVATE
    CGAL::CGAL CGAL::CGAL_Qt5
    Qt5::OpenGL Qt5::Gui)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Circular_kernel_3)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(Circular_kernel_3)

else()

  message(STATUS "NOTICE: This demo requires CGAL, and Qt5, and will not be compiled.")

endif()
