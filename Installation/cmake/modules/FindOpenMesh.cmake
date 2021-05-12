#This modules tries to find OpenMesh
# Once done this will define
#
#  OpenMesh_FOUND - system has OpenMesh
#  OPENMESH_INCLUDE_DIR - OpenMesh include directory
#  OPENMESH_LIBRARIES - OpenMesh libraries
#

find_package(OpenMesh NO_MODULE QUIET)

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

find_package_handle_standard_args(OpenMesh
  REQUIRED_VARS OPENMESH_INCLUDE_DIR OPENMESH_LIBRARIES
  FOUND_VAR OpenMesh_FOUND
  )

if(OpenMesh_FOUND AND NOT TARGET OpenMesh::OpenMesh)
  add_library(OpenMesh::OpenMesh UNKNOWN IMPORTED)

  if(TARGET OpenMeshCore)
    target_link_libraries(OpenMesh::OpenMesh PUBLIC OpenMeshCore)
    return()
  endif()

  set_target_properties(OpenMesh::OpenMesh PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "CGAL_USE_OPENMESH;NOMINMAX;_USE_MATH_DEFINES"
    INTERFACE_INCLUDE_DIRECTORIES  "${OPENMESH_INCLUDE_DIR}")

  if(OPENMESH_LIBRARY_RELEASE)
    set_property(TARGET OpenMesh::OpenMesh APPEND PROPERTY
      IMPORTED_CONFIGURATIONS RELEASE)
    set_target_properties(OpenMesh::OpenMesh PROPERTIES
      IMPORTED_LOCATION_RELEASE "${OPENMESH_LIBRARY_RELEASE}")
  endif()

  if(OPENMESH_LIBRARY_DEBUG)
    set_property(TARGET OpenMesh::OpenMesh APPEND PROPERTY
      IMPORTED_CONFIGURATIONS DEBUG)
    set_target_properties(OpenMesh::OpenMesh PROPERTIES
      IMPORTED_LOCATION_DEBUG "${OPENMESH_LIBRARY_DEBUG}")
  endif()
endif()
