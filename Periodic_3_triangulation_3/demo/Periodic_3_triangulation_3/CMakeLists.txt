
# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Periodic_3_triangulation_3_demo )

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

# Find CGAL
find_package(CGAL COMPONENTS Qt5)
include( ${CGAL_USE_FILE} )

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Xml OpenGL Help Core)

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
if(Qt5_FOUND)
  add_definitions(-DQT_NO_KEYWORDS)
  find_package(QGLViewer)
endif(Qt5_FOUND)

if (CGAL_FOUND AND CGAL_Qt5_FOUND AND OPENGL_FOUND AND Qt5_FOUND AND QGLVIEWER_FOUND AND TARGET Qt5::qcollectiongenerator)

   include_directories( ${QGLVIEWER_INCLUDE_DIR} )

  # UI files (Qt Designer files)
  qt5_wrap_ui ( UI_FILES MainWindow.ui )

  # qrc files (resource files)
  qt5_add_resources( RESOURCE_FILES ./Periodic_3_triangulation_3.qrc )

  # use the Qt MOC preprocessor on classes that derive from QObject
  qt5_generate_moc( "Scene.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_Scene.cpp" )
  qt5_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_MainWindow.cpp" )
  qt5_generate_moc( "Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_Viewer.cpp" )

  if(DEFINED QT_QCOLLECTIONGENERATOR_EXECUTABLE)
  else()
    set(QT_QCOLLECTIONGENERATOR_EXECUTABLE qcollectiongenerator)
  endif()

  # generate QtAssistant collection file
  add_custom_command ( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Periodic_3_triangulation_3.qhc
      DEPENDS Periodic_3_triangulation_3.qhp Periodic_3_triangulation_3.qhcp
      COMMAND Qt5::qcollectiongenerator
                 ${CMAKE_CURRENT_SOURCE_DIR}/Periodic_3_triangulation_3.qhcp
                 -o ${CMAKE_CURRENT_BINARY_DIR}/Periodic_3_triangulation_3.qhc
  )

  # Make sure the compiler can find generated .moc files
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories( ${QT_INCLUDE_DIR} )
  include_directories( ../../include .)

  # The executable itself
  add_executable  ( periodic_3_triangulation_3_demo
                    Scene.cpp moc_Scene.cpp
                    Viewer.cpp moc_Viewer.cpp
                    periodic_3_triangulation_3_demo.cpp
                    MainWindow.ui moc_MainWindow.cpp
                    ${UI_FILES} ${RESOURCE_FILES} Periodic_3_triangulation_3.qhc)

  qt5_use_modules(periodic_3_triangulation_3_demo Xml Help OpenGL)

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS periodic_3_triangulation_3_demo)

  # Link the executable to CGAL and third-party libraries
  target_link_libraries(periodic_3_triangulation_3_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else ()

  set(PERIODIC_TRIANGULATION_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "the CGAL Qt5 library, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "Qt5, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "OpenGL, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "QGLViewer, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if (NOT TARGET Qt5::qcollectiongenerator)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "qcollectiongenerator, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${PERIODIC_TRIANGULATION_MISSING_DEPS}and will not be compiled.")

endif ()
