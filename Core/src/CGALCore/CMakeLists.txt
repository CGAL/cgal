project (CGAL_Core)
   
# Minimal version of CMake:
cmake_minimum_required(VERSION 2.6.2)

# Tested version:
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

if ( NOT CGAL_FOUND )
  find_package(CGAL REQUIRED)
endif()

include(${CGAL_USE_FILE})

include_directories( ${CGAL_3RD_PARTY_INCLUDE_DIRS} )

build_cgal_library(CGAL_Core CGALCore "")

add_dependencies( CGAL_Core CGAL )

add_definitions ( ${CGAL_3RD_PARTY_DEFINITIONS} )
  
target_link_libraries( CGAL_Core ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} )
