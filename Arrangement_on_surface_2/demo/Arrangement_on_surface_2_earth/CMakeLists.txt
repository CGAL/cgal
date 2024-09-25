# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.29)
project(Arrangement_on_surface_2_earth_Demo)

if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(Qt6 QUIET COMPONENTS Core Gui OpenGL OpenGLWidgets Widgets Xml)
find_package(CGAL COMPONENTS Qt6)
find_package(nlohmann_json QUIET 3.9)

if (NOT CGAL_FOUND OR NOT CGAL_Qt6_FOUND OR NOT Qt6_FOUND OR NOT Boost_FOUND OR NOT nlohmann_json_FOUND)
  if (NOT CGAL_FOUND)
    set(MISSING_DEPS "the CGAL library, ${MISSING_DEPS}")
  endif()
  if (NOT CGAL_Qt6_FOUND)
    set(MISSING_DEPS "the CGAL Qt6 component, ${MISSING_DEPS}")
  endif()
  if (NOT Qt6_FOUND)
    set(MISSING_DEPS "the Qt6 library, ${MISSING_DEPS}")
  endif()
  if (NOT Boost_FOUND)
    set(MISSING_DEPS "the Boost library, ${MISSING_DEPS}")
  endif()
  if (NOT nlohmann_json_FOUND)
    set(MISSING_DEPS "JSON for Modern C++ 3.9+ (know as nlohmann_json), ${MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This project requires ${MISSING_DEPS} and will not be compiled.")
  return()
endif()



add_definitions(-DQT_NO_VERSION_TAGGING)

# AOS
file(GLOB source_files_aos Aos.h Aos.cpp Aos_defs.h
  Aos_triangulator.h Aos_triangulator.cpp)
source_group("Aos" FILES ${source_files_aos})

# GIS
file(GLOB source_files_gis Kml_reader.h Kml_reader.cpp)
source_group("GIS" FILES ${source_files_gis})

# GRAPHICS
file(GLOB source_files_graphics Camera.h Camera.cpp
  Shader_program.h Shader_program.cpp)
source_group("Graphics" FILES ${source_files_graphics})

# GRAPHICS-GEOMETRY (Graphics-related)
file(GLOB source_files_graphics_geometry
  Line_strips.h Line_strips.cpp
  Single_vertex.h Single_vertex.cpp
  Sphere.h Sphere.cpp
  Triangles.h Triangles.cpp
  Vertices.h Vertices.cpp
  World_coordinate_axes.h World_coordinate_axes.cpp)
source_group("Graphics_Geometry" FILES ${source_files_graphics_geometry})

# GUI
file(GLOB source_files_gui
  Camera_manip.h Camera_manip.cpp
  Camera_manip_rot.h Camera_manip_rot.cpp
  Camera_manip_rot_bpa.h Camera_manip_rot_bpa.cpp
  Camera_manip_zoom.h Camera_manip_zoom.cpp
  GUI_country_pick_handler.h GUI_country_pick_handler.cpp
  GUI_event_handler.h GUI_event_handler.cpp)
source_group("GUI" FILES ${source_files_gui})

#SOURCE FILES (NOT CATEGORIZED YET)
file(GLOB source_files Common_defs.h earth.cpp Timer.h
  Main_widget.h Main_widget.cpp
  Message_manager.h Message_manager.cpp
  Tools.h Tools.cpp
  Verification.h Verification.cpp)
source_group("Source Files" FILES ${source_files})

set(CMAKE_AUTOMOC ON)

qt_add_executable(earth
  ${source_files_aos}
  ${source_files_gis}
  ${source_files_graphics}
  ${source_files_graphics_geometry}
  ${source_files_gui}
  ${source_files}
)

set_target_properties(earth PROPERTIES
  WIN32_EXECUTABLE TRUE
  MACOSX_BUNDLE TRUE
)

if (MSVC)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:CONSOLE")
endif (MSVC)

target_link_libraries(earth PRIVATE
  Qt6::Core
  Qt6::Gui
  Qt6::OpenGL
  Qt6::OpenGLWidgets
  Qt6::Widgets
  Qt6::Xml
  CGAL::CGAL
  CGAL::Data
  nlohmann_json::nlohmann_json)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# install(TARGETS earth
#   RUNTIME DESTINATION "${INSTALL_EXAMPLEDIR}"
#   BUNDLE DESTINATION "${INSTALL_EXAMPLEDIR}"
#   LIBRARY DESTINATION "${INSTALL_EXAMPLEDIR}"
# )

add_to_cached_list(CGAL_EXECUTABLE_TARGETS earth)
include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
cgal_add_compilation_test(earth)
