# Try to find the CGAL libraries
# CGAL_FOUND             - system has CGAL lib
# CGAL_INCLUDE_DIR       - the CGAL include directory 
# CGAL_LIBRARIES         - the CGAL libraries 
# CGAL_LIBRARIESDIR      - the CGAL libraries directory
# CGAL_CMAKE_MODULE_PATH - the directory containing CGAL related CMake stuff

# CGAL_ROOT              - set to the enviroment variable CGALROOT, IF defined

include(FindPackageHandleStandardArgs)
include(GeneratorSpecificSettings)

if (CGAL_INCLUDE_DIRS AND CGAL_LIBRARIES )
  # Already in cache, be silent
  set(CGAL_FIND_QUIETLY TRUE)
endif()

set(CGAL_ROOT $ENV{CGALROOT})

if(CGAL_ROOT)
  FILE(TO_CMAKE_PATH ${CGAL_ROOT} CGAL_ROOT)
endif()

find_path(CGAL_INCLUDE_DIR 
          NAMES CGAL/basic.h
          PATHS ../../include ${CGAL_ROOT}/include
      	  DOC "The directories containing include files for CGAL"
         )
         
find_path(CGAL_CMAKE_MODULE_PATH CGALcommon.cmake
          PATHS ../../cmake/modules ${CGAL_ROOT}/cmake/modules ${CMAKE_MODULE_PATH}
      	  DOC "The directories containing CGAL CMake stuff"
         )

if ( "${CMAKE_GENERATOR}" MATCHES "Visual Studio" )

  find_path(CGAL_LIBRARIES_DIR 
            NAMES "cgal${TOOLSET}-mt.lib" "cgal${TOOLSET}-mt-gd.lib" "cgal${TOOLSET}-mt-o.lib" "cgal${TOOLSET}-mt-g.lib"
            PATHS ../../lib ${CGAL_ROOT}/lib
       	    DOC "Path to CGAL and third-party libraries"
           ) 

  set(CGAL_LIBRARIES "") 
  
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(CGAL "DEFAULT_MSG" CGAL_LIBRARIES_DIR CGAL_INCLUDE_DIR CGAL_CMAKE_MODULE_PATH)
  
else()

  find_library(CGAL_LIBRARIES NAMES CGAL
               PATHS ../../lib ${CGAL_ROOT}/lib
         	     DOC "CGAL libraries"
              )
              
  get_filename_component(CGAL_LIBRARIES_DIR ${CGAL_LIBRARIES} PATH)
            
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(CGAL "DEFAULT_MSG" CGAL_LIBRARIES CGAL_INCLUDE_DIR CGAL_CMAKE_MODULE_PATH)
endif()
         
