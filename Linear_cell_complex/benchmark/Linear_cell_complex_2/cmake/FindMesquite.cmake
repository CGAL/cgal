# Find Mesquite. If found, this will define
#
# MESQUITE_FOUND       - Successfully found Mesquite
# MESQUITE_INCLUDE_DIR - Mesquite include directory
# MESQUITE_LIBRARY_DIR - Mesquite library directory
# MESQUITE_LIBRARIES   - Mesquite libraries
#

if(DEFINED MESQUITE_INCLUDE_DIR)
  set(MESQUITE_FIND_QUIETLY TRUE)
else()

  find_path(MESQUITE_INCLUDE_DIR Mesquite_all_headers.hpp
    PATHS
    /usr/local/include
    /usr/include
    $ENV{MESQUITE_DIR}/include
    ${MESQUITE_DIR}/include
    )

  if(NOT ${MESQUITE_INCLUDE_DIR} STREQUAL "MESQUITE_INCLUDE_DIR-NOTFOUND")

    message(STATUS "Found Mesquite: " ${MESQUITE_INCLUDE_DIR})
    set(MESQUITE_FOUND true)

    set(MESQUITE_LIBRARY_DIR "${MESQUITE_INCLUDE_DIR}/../lib/" CACHE PATH "Mesquite library directory")
    set(MESQUITE_LIBRARIES "mesquite" CACHE STRING "Mesquite libraries")

  else()
    set(MESQUITE_FOUND FALSE)
  endif()
endif()
