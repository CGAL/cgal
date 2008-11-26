# - Try to find QGLViewer
# Once done this will define
#
#  QGLVIEWER_FOUND - system has QGLViewer
#  QGLVIEWER_INCLUDE_DIR - the QGLViewer include directory
#  QGLVIEWER_LIBRARY - Link these to use QGLViewer
#  QGLVIEWER_DEFINITIONS - Compiler switches required for using QGLViewer
#

find_path(QGLVIEWER_INCLUDE_DIR 
          NAMES QGLViewer/qglviewer.h
          PATHS /usr/include
                /usr/local/include
                ENV QGLVIEWERROOT 
         )

find_library(QGLVIEWER_LIBRARY_RELEASE 
             NAMES QGLViewer QGLViewer2
             PATHS /usr/lib
                   /usr/local/lib
                   ENV QGLVIEWERROOT
                   ENV LD_LIBRARY_PATH
                   ENV LIBRARY_PATH
             PATH_SUFFIXES QGLViewer QGLViewer/release
            )

find_library(QGLVIEWER_LIBRARY_DEBUG
             NAMES dQGLViewer dQGLViewer2
             PATHS /usr/lib
                   /usr/local/lib
                   ENV QGLVIEWERROOT
                   ENV LD_LIBRARY_PATH
                   ENV LIBRARY_PATH
             PATH_SUFFIXES QGLViewer QGLViewer/debug      
            )

set( QGLVIEWER_LIBRARIES optimized ${QGLVIEWER_LIBRARY_RELEASE} debug ${QGLVIEWER_LIBRARY_DEBUG} )

set(QGLVIEWER_LIBRARY ${QGLVIEWER_LIBRARIES} CACHE FILEPATH "The QGLViewer library")

IF(QGLVIEWER_INCLUDE_DIR AND QGLVIEWER_LIBRARY)
   SET(QGLVIEWER_FOUND TRUE)
ENDIF(QGLVIEWER_INCLUDE_DIR AND QGLVIEWER_LIBRARY)

IF(QGLVIEWER_FOUND)
  IF(NOT QGLViewer_FIND_QUIETLY)
    MESSAGE(STATUS "Found QGLViewer: ${QGLVIEWER_LIBRARY}")
  ENDIF(NOT QGLViewer_FIND_QUIETLY)
ELSE(QGLVIEWER_FOUND)
  IF(QGLViewer_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find QGLViewer")
  ENDIF(QGLViewer_FIND_REQUIRED)
ENDIF(QGLVIEWER_FOUND)

