# This file must be included in your CMakeLists.txt to use the basic viewer
# You need to link the libraries in your executable by using
# TARGET_LINK_LIBRARIES( myexec ${BASIC_VIEWER_LIBRARIES})

if ( NOT CGAL_FOUND OR NOT CGAL_Qt5_FOUND)
  message(STATUS "NOTICE: Libraries for basic viewer not found "
                 "(CGAL, Qt5, QGLViewer).")
endif( NOT CGAL_FOUND OR NOT CGAL_Qt5_FOUND)

include( ${CGAL_USE_FILE} )
set(CMAKE_INCLUDE_CURRENT_DIR ON)

FIND_PACKAGE(Qt5        REQUIRED COMPONENTS OpenGL Xml)
find_package(QGLViewer  REQUIRED)
find_package(OpenGL     REQUIRED)

add_definitions(${QT_DEFINITIONS})
add_definitions(-DQT_NO_KEYWORDS)

include_directories( ${QGLVIEWER_INCLUDE_DIR} )
add_definitions(${QGLVIEWER_DEFINITIONS})

set (BASIC_VIEWER_LIBRARIES ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES}
                          ${OPENGL_gl_LIBRARY}) # ${OPENGL_glu_LIBRARY}

set(BASIC_VIEWER_MODULES Xml OpenGL)

ADD_DEFINITIONS("-DCGAL_USE_BASIC_VIEWER")
message(STATUS "Libraries for lcc_viewer found. You need to link them "
               "in your executable by using "
               "TARGET_LINK_LIBRARIES( myexec \${BASIC_VIEWER_LIBRARIES})")

set(USE_BASIC_VIEWER true)
