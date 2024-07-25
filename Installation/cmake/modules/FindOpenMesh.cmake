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
            PATH_SUFFIXES src include
            DOC "The directory containing the OpenMesh header files WITHOUT the OpenMesh prefix"
           )

  find_library(OPENMESH_CORE_LIBRARY_RELEASE NAMES "OpenMeshCore"
               HINTS ENV OPENMESH_LIB_DIR
                     ENV OPENMESH_DIR
               PATH_SUFFIXES lib
               DOC "Path to the OpenMeshCore library"
              )

  find_library(OPENMESH_CORE_LIBRARY_DEBUG NAMES "OpenMeshCored"
               HINTS ENV OPENMESH_LIB_DIR
                     ENV OPENMESH_DIR
               PATH_SUFFIXES lib
               DOC "Path to the OpenMeshCored library"
              )

  find_library(OPENMESH_TOOLS_LIBRARY_RELEASE NAMES "OpenMeshTools"
               HINTS ENV OPENMESH_LIB_DIR
                     ENV OPENMESH_DIR
               PATH_SUFFIXES lib
               DOC "Path to the OpenMeshTools library"
              )

  find_library(OPENMESH_TOOLS_LIBRARY_DEBUG NAMES "OpenMeshToolsd"
               HINTS ENV OPENMESH_LIB_DIR
                     ENV OPENMESH_DIR
               PATH_SUFFIXES lib
               DOC "Path to the OpenMeshToolsd library"
              )

  #select configuration depending on platform (optimized... on windows)
  include(SelectLibraryConfigurations)
  select_library_configurations( OPENMESH_TOOLS )
  select_library_configurations( OPENMESH_CORE )

  set(OPENMESH_LIBRARIES ${OPENMESH_CORE_LIBRARY} ${OPENMESH_TOOLS_LIBRARY} )
  set(OPENMESH_INCLUDE_DIRS ${OPENMESH_INCLUDE_DIR} )

  include( FindPackageHandleStandardArgs )
  find_package_handle_standard_args(OpenMesh
    REQUIRED_VARS OPENMESH_INCLUDE_DIRS OPENMESH_LIBRARIES
    FOUND_VAR OpenMesh_FOUND
  )

  #target OpenMesh::OpenMesh
  if(OpenMesh_FOUND AND NOT TARGET OpenMesh::OpenMesh)
    add_library(CGAL_OpenMesh::CGAL_OpenMesh UNKNOWN IMPORTED)

    if(TARGET OpenMeshCore)
      target_link_libraries(CGAL_OpenMesh::CGAL_OpenMesh INTERFACE OpenMeshCore)
      return()
    endif()

    set_target_properties(CGAL_OpenMesh::CGAL_OpenMesh PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "CGAL_USE_OPENMESH;NOMINMAX;_USE_MATH_DEFINES"
      INTERFACE_INCLUDE_DIRECTORIES  "${OPENMESH_INCLUDE_DIRS}")

    if(OPENMESH_CORE_LIBRARY_RELEASE)
      set_property(TARGET CGAL_OpenMesh::CGAL_OpenMesh APPEND PROPERTY
        IMPORTED_CONFIGURATIONS RELEASE)
      set_target_properties(CGAL_OpenMesh::CGAL_OpenMesh PROPERTIES
        IMPORTED_LOCATION_RELEASE "${OPENMESH_CORE_LIBRARY_RELEASE}")
    endif()

    if(OPENMESH_CORE_LIBRARY_DEBUG)
      set_property(TARGET CGAL_OpenMesh::CGAL_OpenMesh APPEND PROPERTY
        IMPORTED_CONFIGURATIONS DEBUG)
      set_target_properties(CGAL_OpenMesh::CGAL_OpenMesh PROPERTIES
        IMPORTED_LOCATION_DEBUG "${OPENMESH_CORE_LIBRARY_DEBUG}")
    endif()
  endif()

endif()

if(OpenMesh_FOUND)
  message(STATUS "OpenMesh found")
endif()