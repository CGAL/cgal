# This is the CMake script for compiling a CGAL application.

project (GraphicsView)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

find_package(CGAL REQUIRED COMPONENTS Qt5)

include(${CGAL_USE_FILE})

find_package(Qt5 COMPONENTS Xml Script OpenGL)

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )

  add_definitions(-DQT_NO_KEYWORDS)

  add_executable  ( min min.cpp )

  qt5_use_modules(min Xml Script OpenGL)

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS min )

  # Link with Qt libraries
  target_link_libraries( min ${QT_LIBRARIES} )
  # Link with CGAL
  target_link_libraries( min ${CGAL_LIBRARIES}  ${CGAL_3RD_PARTY_LIBRARIES})

else()

  message(STATUS "NOTICE: This demo requires CGAL and Qt5, and will not be compiled.")

endif()
