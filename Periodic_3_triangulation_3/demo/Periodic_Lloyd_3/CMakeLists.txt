# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Periodic_Lloyd_3)

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
  cmake_policy(VERSION 2.8.4)
else()
  cmake_policy(VERSION 2.6)
endif()

set( QT_USE_QTXML    TRUE )
set( QT_USE_QTMAIN   TRUE )
set( QT_USE_QTSCRIPT  TRUE )
set( QT_USE_QTOPENGL  TRUE )

#New for Qt5 version !
option (USE_QT5 "Use Qt5 version instead of Qt4" ON)

if(USE_QT5)

	find_package(CGAL COMPONENTS Qt5)

	include(${CGAL_USE_FILE})
	
	find_package(Qt5)

else (USE_QT5)

	find_package(CGAL COMPONENTS Qt4)
	include(${CGAL_USE_FILE})

	find_package(Qt4)

	if(CGAL_FOUND AND QT4_FOUND)
		include(${QT_USE_FILE})
		message("Qt4 found")
		set(QT4 TRUE)
	endif()

endif (USE_QT5)

include(QtChoice)

find_package(OpenGL)

if(QT4 OR QT5)
  find_package(QGLViewer)
endif(QT4 OR QT5)

if ( CGAL_FOUND AND (QT4 OR QT5) AND OPENGL_FOUND AND QGLVIEWER_FOUND AND QT_QCOLLECTIONGENERATOR_EXECUTABLE )


  include_directories (${QGLVIEWER_INCLUDE_DIR})
  include_directories (BEFORE ../../include ./ )

  # ui file, created wih Qt Designer
  qt_wrap_ui( uis MainWindow.ui )

  # qrc files (resources files, that contain icons, at least)
  qt_add_resources ( RESOURCE_FILES ./Periodic_Lloyd_3.qrc )

  qt_automoc( MainWindow.cpp Viewer.cpp)

  if(DEFINED QT_QCOLLECTIONGENERATOR_EXECUTABLE)
  else()
    set(QT_QCOLLECTIONGENERATOR_EXECUTABLE qcollectiongenerator)
  endif()

  # generate QtAssistant collection file
  add_custom_command ( OUTPUT Periodic_Lloyd_3.qhc
      DEPENDS Periodic_Lloyd_3.qhp Periodic_Lloyd_3.qhcp
      COMMAND ${QT_QCOLLECTIONGENERATOR_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/Periodic_Lloyd_3.qhcp
                -o ${CMAKE_CURRENT_BINARY_DIR}/Periodic_Lloyd_3.qhc
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_executable  ( Periodic_Lloyd_3 
                    Periodic_Lloyd_3.qhc 
                    Periodic_Lloyd_3.cpp 
                    MainWindow.cpp 
                    Viewer.cpp 
                    ${uis}  
                    ${RESOURCE_FILES} )

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Periodic_Lloyd_3 )

  target_link_libraries( Periodic_Lloyd_3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  target_link_libraries( Periodic_Lloyd_3 ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES} )
  target_link_libraries( Periodic_Lloyd_3 ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} )

else ( CGAL_FOUND AND (QT4 OR QT5) AND OPENGL_FOUND AND QGLVIEWER_FOUND AND QT_QCOLLECTIONGENERATOR_EXECUTABLE )

  set(PERIODIC_LLOYD_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "the CGAL Qt library, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()      	   

  if(NOT QT4 AND NOT USE_QT5)
    set(PERIODIC_LLOYD_MISSING_DEPS "Qt4, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()  

  if(NOT QT5 AND USE_QT5)
    set(PERIODIC_LLOYD_MISSING_DEPS "Qt5, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()         	   

  if(NOT OPENGL_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "OpenGL, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()      	   

  if(NOT QGLVIEWER_FOUND)
    set(PERIODIC_LLOYD_MISSING_DEPS "QGLViewer, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()      	   

  if(NOT QT_QCOLLECTIONGENERATOR_EXECUTABLE)
    set(PERIODIC_LLOYD_MISSING_DEPS "qcollectiongenerator, ${PERIODIC_LLOYD_MISSING_DEPS}")
  endif()


  message(STATUS "NOTICE: This demo requires ${PERIODIC_LLOYD_MISSING_DEPS}and will not be compiled.")

endif ( CGAL_FOUND AND (QT4 OR QT5) AND OPENGL_FOUND AND QGLVIEWER_FOUND AND QT_QCOLLECTIONGENERATOR_EXECUTABLE )
