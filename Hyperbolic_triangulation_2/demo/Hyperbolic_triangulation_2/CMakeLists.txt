cmake_minimum_required(VERSION 3.1...3.23)
project(Hyperbolic_triangulation_2_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Core Qt6)
find_package(LEDA QUIET)

# Find Qt6 itself
find_package(Qt6 QUIET COMPONENTS Widgets OpenGLWidgets)

if(CGAL_Qt6_FOUND
   AND Qt6_FOUND
   AND (CGAL_Core_FOUND OR LEDA_FOUND))
  # ui files, created with Qt Designer
  qt6_wrap_ui(UIS HDT2.ui)

  qt6_add_resources(RESOURCE_FILES resources/Delaunay_triangulation_2.qrc)

  # cpp files
  add_executable ( HDT2 HDT2.cpp ${CGAL_Qt6_RESOURCE_FILES} ${RESOURCE_FILES} ${UIS})
  target_include_directories(HDT2 PRIVATE ./ ./include)
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS HDT2 )
  target_link_libraries ( HDT2 CGAL::CGAL CGAL::CGAL_Qt6 Qt6::Widgets)
  if(CGAL_Core_FOUND)
    target_link_libraries ( HDT2 CGAL::CGAL_Core)
  else()
    target_link_libraries ( HDT2 ${LEDA_LIBRARIES})
  endif()

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test( HDT2 )
else()
  message("NOTICE: This demo requires CGAL_Core (or LEDA), and Qt6 and will not be compiled.")
endif()
