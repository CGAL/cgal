# Copyright (C) 2022 The Qt Company Ltd.
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required(VERSION 3.16)
project(earth LANGUAGES CXX)

if(NOT DEFINED INSTALL_EXAMPLESDIR)
    set(INSTALL_EXAMPLESDIR "examples")
endif()

set(INSTALL_EXAMPLEDIR "${INSTALL_EXAMPLESDIR}/opengl/earth")

find_package(Qt6 REQUIRED COMPONENTS Core Gui OpenGL OpenGLWidgets Widgets)

qt_standard_project_setup()

qt_add_executable(earth
    main.cpp
    mainwidget.h      mainwidget.cpp
    Camera.h          Camera.cpp
    Common_defs.h
    Shader_program.h  Shader_program.cpp
    Sphere.h          Sphere.cpp
    Tools.h
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
)


# Resources:
#set(shaders_resource_files
#    "vertex_shader.glsl"
#)
#
#qt6_add_resources(earth "shaders"
#    PREFIX
#        "/"
#    FILES
#        ${shaders_resource_files}
#)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


install(TARGETS earth
    RUNTIME DESTINATION "${INSTALL_EXAMPLEDIR}"
    BUNDLE DESTINATION "${INSTALL_EXAMPLEDIR}"
    LIBRARY DESTINATION "${INSTALL_EXAMPLEDIR}"
)
