# This is the CMake script for compiling the AABB tree demo.

cmake_minimum_required(VERSION 3.1...3.23)
project(AABB_tree_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Find CGAL and CGAL Qt6
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

# Find Qt6 itself
find_package(Qt6 QUIET COMPONENTS Gui OpenGL)

if(CGAL_Qt6_FOUND AND Qt6_FOUND)

  add_definitions(-DQT_NO_KEYWORDS)

  # Instruct CMake to run moc/ui/rcc automatically when needed.
  set(CMAKE_AUTOMOC ON)
  set(CMAKE_AUTOUIC ON)
  set(CMAKE_AUTORCC ON)

  qt_add_executable(
    AABB_demo AABB_demo.cpp Scene.cpp benchmarks.cpp Viewer.cpp MainWindow.cpp
    MainWindow.ui AABB_demo.qrc
  )
  # Link with Qt libraries
  target_link_libraries(AABB_demo PRIVATE Qt6::Gui Qt6::OpenGL
                                          CGAL::CGAL CGAL::CGAL_Qt6)

  add_to_cached_list(CGAL_EXECUTABLE_TARGETS AABB_demo)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(AABB_demo)

else(CGAL_Qt6_FOUND AND Qt6_FOUND)

  set(AABB_MISSING_DEPS "")

  if(NOT CGAL_Qt6_FOUND)
    set(AABB_MISSING_DEPS "CGAL_Qt6, ${AABB_MISSING_DEPS}")
  endif()

  if(NOT Qt6_FOUND)
    set(AABB_MISSING_DEPS "Qt6, ${AABB_MISSING_DEPS}")
  endif()

  message("NOTICE: This demo requires ${AABB_MISSING_DEPS}, and will not be compiled.")

endif(CGAL_Qt6_FOUND AND Qt6_FOUND)
