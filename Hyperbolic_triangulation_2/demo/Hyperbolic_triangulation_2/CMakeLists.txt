project (Hyperbolic_triangulation_2_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON) 

cmake_minimum_required(VERSION 3.1...3.13)

find_package(CGAL QUIET COMPONENTS Core Qt5)

find_package(Qt5 QUIET COMPONENTS Widgets)

include_directories (BEFORE include)

if(CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND CGAL_Core_FOUND)
  # ui files, created with Qt Designer
  qt5_wrap_ui( UIS HDT2.ui )
  
  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources ( RESOURCE_FILES resources/File.qrc resources/Input.qrc resources/Delaunay_triangulation_2.qrc )
  
  # cpp files
  add_executable ( HDT2_demo HDT2.cpp ${RESOURCE_FILES} ${UIS})
  target_link_libraries ( HDT2_demo CGAL::CGAL CGAL::CGAL_Qt5 CGAL::CGAL_Core Qt5::Widgets)
else()
  message(STATUS "NOTICE: This demo requires CGAL, Boost and Qt5 and will not be compiled.")
endif()
