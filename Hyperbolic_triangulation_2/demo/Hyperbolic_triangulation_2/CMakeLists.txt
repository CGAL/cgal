cmake_minimum_required(VERSION 3.1...3.15)
project(Hyperbolic_triangulation_2_Demo)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(CGAL REQUIRED QUIET OPTIONAL_COMPONENTS Core Qt5)
find_package(LEDA QUIET)

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS OpenGL Gui)

if(CGAL_Qt5_FOUND
   AND Qt5_FOUND
   AND (CGAL_Core_FOUND OR LEDA_FOUND))
  include_directories(BEFORE ./ ./include)
  # ui files, created with Qt Designer
  qt5_wrap_ui(UIS HDT2.ui)

  qt5_add_resources(RESOURCE_FILES resources/Delaunay_triangulation_2.qrc)

  # cpp files
  add_executable(HDT2 HDT2.cpp ${CGAL_Qt5_RESOURCE_FILES} ${RESOURCE_FILES}
                      ${UIS})
  add_to_cached_list(CGAL_EXECUTABLE_TARGETS HDT2)
  target_link_libraries(HDT2 CGAL::CGAL CGAL::CGAL_Qt5 CGAL::CGAL_Core
                        Qt5::Widgets)
else()
  message(
    STATUS
      "NOTICE: This demo requires CGAL, CGAL_Core (or LEDA), and Qt5 and will not be compiled."
  )
endif()
