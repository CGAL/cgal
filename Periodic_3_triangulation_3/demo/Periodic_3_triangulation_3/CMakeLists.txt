# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Periodic_3_triangulation_3_demo ) 

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

# Find CGAL
find_package(CGAL)
include( ${CGAL_USE_FILE} )

# Find Qt4 itself
set ( QT_USE_QTXML    TRUE )
set ( QT_USE_QTMAIN   TRUE )
set ( QT_USE_QTOPENGL TRUE )
find_package(Qt4)

# Find OpenGL
find_package(OpenGL)

# Find QGLViewer
if(QT4_FOUND)
  include( ${QT_USE_FILE} )
  find_package(QGLViewer)
endif(QT4_FOUND)

if (CGAL_FOUND AND OPENGL_FOUND AND QT4_FOUND AND QGLVIEWER_FOUND) 

   include_directories( ${QGLVIEWER_INCLUDE_DIR} )

  # UI files (Qt Designer files)
  qt4_wrap_ui ( UI_FILES MainWindow.ui )

  # qrc files (resource files)
  qt4_add_resources( RESOURCE_FILES ./Periodic_3_triangulation_3.qrc )

  # use the Qt MOC preprocessor on classes that derive from QObject
  qt4_generate_moc( "Scene.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_Scene.cpp" )
  qt4_generate_moc( "MainWindow.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_MainWindow.cpp" )

  if(DEFINED QT_QCOLLECTIONGENERATOR_EXECUTABLE)
  else()
    set(QT_QCOLLECTIONGENERATOR_EXECUTABLE qcollectiongenerator)
  endif()

  # generate QtAssistant collection file
  add_custom_command ( OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Periodic_3_triangulation_3.qhc
      DEPENDS Periodic_3_triangulation_3.qhp Periodic_3_triangulation_3.qhcp
      COMMAND ${QT_QCOLLECTIONGENERATOR_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/Periodic_3_triangulation_3.qhcp
                 -o ${CMAKE_CURRENT_BINARY_DIR}/Periodic_3_triangulation_3.qhc
  )

  # Make sure the compiler can find generated .moc files
  include_directories(BEFORE ${CMAKE_CURRENT_BINARY_DIR})
  include_directories( ${QT_INCLUDE_DIR} )
  include_directories( ../../include )

  # The executable itself
  add_executable  ( periodic_3_triangulation_3_demo
  		    Scene.cpp moc_Scene.cpp
  		    periodic_3_triangulation_3_demo.cpp
  		    MainWindow.ui moc_MainWindow.cpp 
		    ${UI_FILES} ${RESOURCE_FILES} Periodic_3_triangulation_3.qhc)

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS periodic_3_triangulation_3_demo)

  # Link the executable to CGAL and third-party libraries
  target_link_libraries(periodic_3_triangulation_3_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else (CGAL_FOUND AND OPENGL_FOUND AND QT4_FOUND AND QGLVIEWER_FOUND) 

  set(PERIODIC_TRIANGULATION_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "the CGAL Qt4 library, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT QT4_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "Qt4, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "OpenGL, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  if(NOT QGLVIEWER_FOUND)
    set(PERIODIC_TRIANGULATION_MISSING_DEPS "QGLViewer, ${PERIODIC_TRIANGULATION_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${PERIODIC_TRIANGULATION_MISSING_DEPS}and will not be compiled.")

endif (CGAL_FOUND AND OPENGL_FOUND AND QT4_FOUND AND QGLVIEWER_FOUND) 
