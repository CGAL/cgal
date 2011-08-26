# - Try to find QGLViewer
# Once done this will define
#
#  QGLVIEWER_FOUND - system has QGLViewer
#  QGLVIEWER_INCLUDE_DIR - the QGLViewer include directory
#  QGLVIEWER_LIBRARIES - Link these to use QGLViewer
#  QGLVIEWER_DEFINITIONS - Compiler switches required for using QGLViewer
#

find_path(QGLVIEWER_INCLUDE_DIR 
          NAMES QGLViewer/qglviewer.h
          PATHS /usr/include
                /usr/local/include
                ENV QGLVIEWERROOT 
         )

find_library(QGLVIEWER_LIBRARY_RELEASE 
             NAMES qglviewer-qt4 qglviewer QGLViewer QGLViewer2
             PATHS /usr/lib
                   /usr/local/lib
                   ENV QGLVIEWERROOT
                   ENV LD_LIBRARY_PATH
                   ENV LIBRARY_PATH
             PATH_SUFFIXES QGLViewer QGLViewer/release
            )

find_library(QGLVIEWER_LIBRARY_DEBUG
             NAMES dqglviewer dQGLViewer dQGLViewer2 QGLViewerd2
             PATHS /usr/lib
                   /usr/local/lib
                   ENV QGLVIEWERROOT
                   ENV LD_LIBRARY_PATH
                   ENV LIBRARY_PATH
             PATH_SUFFIXES QGLViewer QGLViewer/debug      
            )

if(QGLVIEWER_LIBRARY_RELEASE)
  if(QGLVIEWER_LIBRARY_DEBUG)
    set(QGLVIEWER_LIBRARIES_ optimized ${QGLVIEWER_LIBRARY_RELEASE} debug ${QGLVIEWER_LIBRARY_DEBUG})
  else()
    set(QGLVIEWER_LIBRARIES_ ${QGLVIEWER_LIBRARY_RELEASE})
  endif()

  set(QGLVIEWER_LIBRARIES ${QGLVIEWER_LIBRARIES_} CACHE FILEPATH "The QGLViewer library")

endif()

IF(QGLVIEWER_INCLUDE_DIR AND QGLVIEWER_LIBRARIES)
   SET(QGLVIEWER_FOUND TRUE)
ENDIF(QGLVIEWER_INCLUDE_DIR AND QGLVIEWER_LIBRARIES)

IF(QGLVIEWER_FOUND)
  IF(NOT QGLViewer_FIND_QUIETLY)
    MESSAGE(STATUS "Found QGLViewer: ${QGLVIEWER_LIBRARIES}")
  ENDIF(NOT QGLViewer_FIND_QUIETLY)
ELSE(QGLVIEWER_FOUND)
  IF(QGLViewer_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find QGLViewer")
  ENDIF(QGLViewer_FIND_REQUIRED)
ENDIF(QGLVIEWER_FOUND)

