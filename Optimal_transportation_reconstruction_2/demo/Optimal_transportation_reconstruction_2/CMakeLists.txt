# This is the CMake script for compiling the Optimal_transportation_reconstruction_2 demo.

project(Otr2_demo)

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
  cmake_policy(SET CMP0043 OLD)
endif()

# Include this package's headers first
include_directories( BEFORE ./ ./include ../../include )

# Find CGAL and CGAL Qt5
find_package(CGAL COMPONENTS Qt5)
include( ${CGAL_USE_FILE} )

# Find Qt5 itself
find_package(Qt5 5.4 QUIET COMPONENTS OpenGL)
if(Qt5_FOUND)
  add_definitions(-DQT_NO_KEYWORDS)
endif(Qt5_FOUND)

# Find OpenGL
find_package(OpenGL)

# Find CImg
find_path(CIMG_INCLUDE_DIR
          NAMES CImg.h
          HINTS ENV CIMG_INC_DIR
          DOC "Path to the header of the CImg library")

if (CIMG_INCLUDE_DIR)
  add_definitions(-DCGAL_USE_CIMG)
  include_directories(${CIMG_INCLUDE_DIR})
  MESSAGE(STATUS "CImg library found, the demo can load point set from image files.")
else()
  MESSAGE(STATUS "CImg library was not found, the demo will not be able to load point set from image files. "
                 "Try setting the environment variable CIMG_INC_DIR to point to the path of the directory containing CImg.h.")
endif()

if(CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND)

  set(
  SRCS 
  glviewer.cpp 
  Otr2_demo.cpp 
  window.cpp 
  render.cpp 
  )

  set( 
  CGAL_Qt5_MOC_FILES 
  moc_dialog_options.cxx 
  moc_glviewer.cxx 
  moc_window.cxx 
  )

  set(
  UIS
  pwsrec.ui 
  options.ui
  )
  
  qt5_wrap_ui( UI_FILES ${UIS} )

  include(AddFileDependencies)

  qt5_generate_moc( "window.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_window.cxx" )
  add_file_dependencies( moc_window.cxx "${CMAKE_CURRENT_SOURCE_DIR}/window.h" )

  qt5_generate_moc( "glviewer.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_glviewer.cxx" )
  add_file_dependencies( moc_glviewer.cxx "${CMAKE_CURRENT_SOURCE_DIR}/glviewer.h" )

  qt5_generate_moc( "dialog_options.h" "${CMAKE_CURRENT_BINARY_DIR}/moc_dialog_options.cxx" )
  add_file_dependencies( moc_dialog_options.cxx "${CMAKE_CURRENT_SOURCE_DIR}/dialog_options.h" )

  qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES pwsrec.qrc )

  add_executable  ( Otr2_demo ${SRCS} ${CGAL_Qt5_MOC_FILES} ${UI_FILES} ${CGAL_Qt5_RESOURCE_FILES})
  
  qt5_use_modules(Otr2_demo OpenGL)

  # Link with Qt libraries
  target_link_libraries( Otr2_demo ${QT_LIBRARIES} )

  # Link with CGAL
  target_link_libraries( Otr2_demo ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

  # Link with OpenGL
  target_link_libraries( Otr2_demo ${OPENGL_gl_LIBRARY} )
    
  add_to_cached_list( CGAL_EXECUTABLE_TARGETS Otr2_demo )

else (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND )

  set(OTR2_MISSING_DEPS "")

  if(NOT CGAL_FOUND)
    set(OTR2_MISSING_DEPS "the CGAL library, ${OTR2_MISSING_DEPS}")
  endif()

  if(NOT CGAL_Qt5_FOUND)
    set(OTR2_MISSING_DEPS "the CGAL Qt5 library, ${OTR2_MISSING_DEPS}")
  endif()

  if(NOT Qt5_FOUND)
    set(OTR2_MISSING_DEPS "Qt5.4, ${OTR2_MISSING_DEPS}")
  endif()

  if(NOT OPENGL_FOUND)
    set(OTR2_MISSING_DEPS "OpenGL, ${OTR2_MISSING_DEPS}")
  endif()

  message(STATUS "NOTICE: This demo requires ${OTR2_MISSING_DEPS}and will not be compiled.")

endif (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND AND OPENGL_FOUND)
