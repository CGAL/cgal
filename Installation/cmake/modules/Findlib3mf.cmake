# - Try to find the LibSSH libraries
# This module defines:
#  LIB3MF_FOUND             - system has lib3mf lib
#  LIB3MF_INCLUDE_DIR       - the lib3mf include directory
#  LIB3MF_LIBRARIES         - Link these to use lib3mf


include(FindPackageHandleStandardArgs)
include(${CMAKE_CURRENT_LIST_DIR}/CGAL_GeneratorSpecificSettings.cmake)

# Is it already configured?
if(NOT LIB3MF_LIBRARIES OR NOT LIB3MF_INCLUDE_DIR)
  find_path(LIB3MF_INCLUDE_DIR NAMES "lib3mf_implicit.hpp"
            HINTS "/usr/include/lib3mf/Bindings/Cpp"
                  ENV LIB3MF_INC_DIR)

  find_library(LIB3MF_LIBRARIES NAMES 3mf
               HINTS "/usr/lib"
                     "usr/lib/x86_64-linux-gnu"
                     "usr/lib64"
                     ENV LIB3MF_LIB_DIR
               PATH_SUFFIXES lib
               DOC "Path to the lib3mf library")
endif()

SET(LIB3MF_FOUND TRUE)
if( NOT LIB3MF_LIBRARIES OR NOT LIB3MF_INCLUDE_DIR)
  SET(LIB3MF_FOUND FALSE)
endif()
