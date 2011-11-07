# This file must be included in your CMakeLists.txt to use cgal_map_viewer_qt.h.
# You need to link the libraries in your executable by using 
# TARGET_LINK_LIBRARIES( myexec ${MAP_VIEWER_LIBRARIES_QT})

find_package(CGAL COMPONENTS Qt4)

if ( CGAL_FOUND AND CGAL_Qt4_FOUND)
  include( ${CGAL_USE_FILE} )

  SET(QT_USE_QTMAIN   TRUE )
  SET(QT_USE_QTOPENGL TRUE) 
  SET(QT_USE_QTXML    TRUE) 
  FIND_PACKAGE(Qt4        REQUIRED)
  find_package(QGLViewer  REQUIRED)
  find_package(OpenGL     REQUIRED)

  INCLUDE(${QT_USE_FILE})
  add_definitions(${QT_DEFINITIONS})
  
  include_directories( ${QGLVIEWER_INCLUDE_DIR} )
  add_definitions(${QGLVIEWER_DEFINITIONS})

  set (MAP_VIEWER_LIBRARIES_QT ${QT_LIBRARIES} ${QGLVIEWER_LIBRARIES}
                               ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
    
  message(STATUS "Libraries for map_viewer found. You need to link them "
                 "in your executable by using "
                 "TARGET_LINK_LIBRARIES( myexec \${MAP_VIEWER_LIBRARIES_QT})")

else(CGAL_FOUND AND CGAL_Qt4_FOUND)
  message(STATUS "NOTICE: Libraries for map_viewer not found "
                 "(CGAL, Qt4, QGLViewer).")
endif(CGAL_FOUND AND CGAL_Qt4_FOUND)
