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

#set( BUILD_SHARED_LIBS ON )

FIND_LIBRARY(QGLVIEWER_LIBRARY 
             NAMES QGLViewer QGLViewer2 dQGLViewer
             PATHS /usr/lib
                   /usr/local/lib
                   ENV QGLVIEWERROOT
                   ENV LD_LIBRARY_PATH
                   ENV LIBRARY_PATH
             PATH_SUFFIXES QGLViewer/release QGLViewer/debug      
            )

IF(NOT QGLVIEWER_LIBRARY)
  FIND_LIBRARY(QGLVIEWER_LIBRARY 
               NAMES 3dviewer
               PATHS /usr/lib
                     /usr/local/lib
                     ENV QGLVIEWERROOT
                     ENV LD_LIBRARY_PATH
                     ENV LIBRARY_PATH
               PATH_SUFFIXES release debug      
              )
ENDIF(NOT QGLVIEWER_LIBRARY)

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

# show the QGLVIEWER_INCLUDE_DIR and QGLVIEWER_LIBRARY variables only in the advanced view
IF(QGLVIEWER_FOUND)
  MARK_AS_ADVANCED(QGLVIEWER_INCLUDE_DIR QGLVIEWER_LIBRARY )
ENDIF(QGLVIEWER_FOUND)

