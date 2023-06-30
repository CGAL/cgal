# Copyright (C) 2022 The Qt Company Ltd.
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required(VERSION 3.16)
project(earth LANGUAGES CXX)

if(NOT DEFINED INSTALL_EXAMPLESDIR)
    set(INSTALL_EXAMPLESDIR "examples")
endif()

set(INSTALL_EXAMPLEDIR "${INSTALL_EXAMPLESDIR}/opengl/earth")

#find_package(shapelib REQUIRED)
find_package(Qt6 REQUIRED COMPONENTS Core Gui OpenGL OpenGLWidgets Widgets Xml)
add_definitions(-DQT_NO_VERSION_TAGGING)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()


set(SHAPELIB_INCLUDE_DIR "" CACHE PATH "DEFINE ME!!!")
#target_include_directories(earth PRIVATE ${SHAPELIB_INCLUDE_DIR})
include_directories(${SHAPELIB_INCLUDE_DIR})
set(SHAPELIB_LIB_DIR "" CACHE PATH "DEFINE ME!!!")
#target_link_directories(earth PRIVATE ${SHAPELIB_LIB_DIR})
link_directories(earth PRIVATE ${SHAPELIB_LIB_DIR})


file(GLOB source_files_geometry   
  Line_strips.h             Line_strips.cpp
  Sphere.h                  Sphere.cpp
  Vertices.h                Vertices.cpp
  World_coordinate_axes.h   World_coordinate_axes.cpp
 )
source_group( "geometry"  FILES ${source_files_geometry} )


qt_standard_project_setup()

qt_add_executable(earth
    Aos.h                     Aos.cpp
    Camera.h                  Camera.cpp
    Common_defs.h
    Geodesic_arcs.h           Geodesic_arcs.cpp
    Kml_reader.h              Kml_reader.cpp
                              main.cpp
    Main_widget.h             Main_widget.cpp
    Shader_program.h          Shader_program.cpp
    Tools.h                   Tools.cpp
    ${source_files_geometry}
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
    shp
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


install(TARGETS earth
    RUNTIME DESTINATION "${INSTALL_EXAMPLEDIR}"
    BUNDLE DESTINATION "${INSTALL_EXAMPLEDIR}"
    LIBRARY DESTINATION "${INSTALL_EXAMPLEDIR}"
)
