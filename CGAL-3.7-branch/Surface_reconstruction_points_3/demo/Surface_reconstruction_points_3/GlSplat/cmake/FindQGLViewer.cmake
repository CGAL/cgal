
if (QGLViewer_INCLUDES AND QGLViewer_LIBRARIES)
  set(QGLViewer_FIND_QUIETLY TRUE)
endif (QGLViewer_INCLUDES AND QGLViewer_LIBRARIES)

find_path(QGLViewer_INCLUDES
  NAMES
  QGLViewer/qglviewer.h
  PATHS
  $ENV{QGLViewerDIR}
  ${INCLUDE_INSTALL_DIR}
)

find_library(QGLViewer_LIBRARIES QGLViewer PATHS $ENV{QGLVIEWER_DIR} ${LIB_INSTALL_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QGLViewer DEFAULT_MSG
                                  QGLViewer_INCLUDES QGLViewer_LIBRARIES)

mark_as_advanced(QGLViewer_INCLUDES QGLViewer_LIBRARIES)
