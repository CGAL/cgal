# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project (Triangulation_3_Demo)

cmake_minimum_required(VERSION 3.1)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()
if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

find_package(CGAL  COMPONENTS Qt5)
include(${CGAL_USE_FILE})


find_package(Qt5 QUIET COMPONENTS OpenGL)

if(Qt5_FOUND)
  add_definitions(-DQT_NO_KEYWORDS)
  find_package(QGLViewer)
endif(Qt5_FOUND)

# Activate concurrency ? (turned OFF by default)
option(CGAL_ACTIVATE_CONCURRENT_TRIANGULATION_3
  "Activate parallelism in Triangulation_3"
  OFF)
  
# And add -DCGAL_CONCURRENT_TRIANGULATION_3 if that option is ON
if( CGAL_ACTIVATE_CONCURRENT_TRIANGULATION_3 )
  add_definitions( -DCGAL_CONCURRENT_TRIANGULATION_3 )
  find_package( TBB REQUIRED )
else( CGAL_ACTIVATE_CONCURRENT_TRIANGULATION_3 )
  option( LINK_WITH_TBB
    "Link with TBB anyway so we can use TBB timers for profiling"
    ON)
  if( LINK_WITH_TBB )
    find_package( TBB )
  endif( LINK_WITH_TBB )
endif()

if ( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND QGLVIEWER_FOUND )

  include_directories (${QGLVIEWER_INCLUDE_DIR})
  include_directories (BEFORE ../../include ./ )

  # ui files, created with Qt Designer
  qt5_wrap_ui( uis MainWindow.ui )

  # qrc files (resources files, that contain icons, at least)
  qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES ./T3_demo.qrc )


  # cpp files
  add_executable  ( T3_demo T3_demo.cpp
                    MainWindow.cpp Viewer.cpp
                    PreferenceDlg.cpp
                    Scene.cpp ${uis} ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES})

  add_to_cached_list( CGAL_EXECUTABLE_TARGETS T3_demo )

  target_link_libraries( T3_demo PRIVATE ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
  target_link_libraries( T3_demo PRIVATE Qt5::OpenGL ${QGLVIEWER_LIBRARIES})
  if(TBB_FOUND)
    CGAL_target_use_TBB(T3_demo)
  endif()

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(T3_demo)

  include(${CGAL_MODULES_DIR}/CGAL_add_test.cmake)
  cgal_add_compilation_test(T3_demo)

else( CGAL_FOUND AND Qt5_FOUND AND QGLVIEWER_FOUND )

  set(TRIANGULATION_3_MISSING_DEPS "")

  if(NOT CGAL_Qt5_FOUND)
    set(TRIANGULATION_3_MISSING_DEPS "the CGAL Qt5 library, ${TRIANGULATION_3_MISSING_DEPS}")
  endif()      	   

  if(NOT Qt5_FOUND)
    set(TRIANGULATION_3_MISSING_DEPS "Qt5, ${TRIANGULATION_3_MISSING_DEPS}")
  endif()      	   

  if(NOT QGLVIEWER_FOUND)
    set(TRIANGULATION_3_MISSING_DEPS "QGLViewer, ${TRIANGULATION_3_MISSING_DEPS}")
  endif()      	   


  message(STATUS "NOTICE: This demo requires ${TRIANGULATION_3_MISSING_DEPS}and will not be compiled.")

endif( CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND QGLVIEWER_FOUND )
