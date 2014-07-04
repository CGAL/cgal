#This modules tries to find OpenMesh
# Once done this will define
#
#  OpenMesh_FOUND - system has OpenMesh
#  OPENMESH_INCLUDE_DIR - OpenMesh include directory
#  OPENMESH_LIBRARIES - OpenMesh libraries
#

# Is it already configured?
if (NOT OpenMesh_FOUND)

 find_path(OPENMESH_INCLUDE_DIR
            NAMES OpenMesh/Core/Mesh/ArrayKernel.hh
            HINTS ENV OPENMESH_INC_DIR
                  ENV OPENMESH_DIR
                  /usr/include
                  /usr/local/include
            PATH_SUFFIXES src
            DOC "The directory containing the OpenMesh header files WITHOUT the OpenMesh prefix"
           )

  find_library(OPENMESH_LIBRARY_RELEASE NAMES "OpenMeshCore"
               HINTS ENV OPENMESH_LIB_DIR
                     ENV OPENMESH_DIR
               PATH_SUFFIXES lib
               DOC "Path to the OpenMeshCore library"
              )

  find_library(OPENMESH_LIBRARY_DEBUG NAMES "OpenMeshCored"
               HINTS ENV OPENMESH_LIB_DIR
                     ENV OPENMESH_DIR
               PATH_SUFFIXES lib
               DOC "Path to the OpenMeshCored library"
              )

  if(OPENMESH_LIBRARY_RELEASE)
    if(OPENMESH_LIBRARY_DEBUG)
      set(OPENMESH_LIBRARIES optimized ${OPENMESH_LIBRARY_RELEASE} debug ${OPENMESH_LIBRARY_DEBUG})
    else()
      set(OPENMESH_LIBRARIES ${OPENMESH_LIBRARY_RELEASE})
    endif()
  endif()
endif()

include( FindPackageHandleStandardArgs )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(OpenMesh
        REQUIRED_VARS OPENMESH_INCLUDE_DIR OPENMESH_LIBRARIES
        FOUND_VAR OpenMesh_FOUND
        )
