# Find OpenMesh. If found, this will define
#
# OPENMESH_FOUND       - Successfully found OpenMesh
# OPENMESH_INCLUDE_DIR - OpenMesh include directory
# OPENMESH_LIBRARIES   - OpenMesh libraries
# OPENMESH_LIBRARY_DIR - OpenMesh library directory
#

if(DEFINED OPENMESH_INCLUDE_DIR)
  set(OPENMESH_FIND_QUIETLY TRUE)
else()

  find_path(OPENMESH_INCLUDE_DIR OpenMesh/Core/Mesh/PolyMeshT.hh
    PATHS
    /usr/local/include
    /usr/include
    $ENV{OPENMESH_DIR}/include
    ${OPENMESH_DIR}/include
    )

  if(DEFINED OPENMESH_INCLUDE_DIR)

    message(STATUS "Found OpenMesh: " ${OPENMESH_INCLUDE_DIR})
    set(OPENMESH_FOUND true)

    if(WIN32)
      set(OPENMESH_LIBRARY_DIR "${OPENMESH_INCLUDE_DIR}/../lib"
        CACHE PATH "OpenMesh library directory")
    else()
      set(OPENMESH_LIBRARY_DIR "${OPENMESH_INCLUDE_DIR}/../lib/OpenMesh"
        CACHE PATH "OpenMesh library directory")
    endif()

    set(OPENMESH_LIBRARIES "OpenMeshCore;OpenMeshTools"
      CACHE STRING "OpenMesh libraries")

  else()
    set(OPENMESH_FOUND FALSE)
  endif()
endif()
