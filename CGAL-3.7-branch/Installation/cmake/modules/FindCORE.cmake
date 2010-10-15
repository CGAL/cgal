# Try to find the CORE libraries
# CORE_FOUND - system has CORE lib
# CORE_INCLUDE_DIR - the CORE include directory
# CORE_LIBRARIES - Libraries needed to use CORE

# TODO: support Windows and MacOSX

# CORE needs GMP
include(CGAL_FindPackageHandleStandardArgs)

if(GMP_FOUND)
  if (CORE_INCLUDE_DIR AND CORE_LIBRARIES)
    # Already in cache, be silent
    set(CORE_FIND_QUIETLY TRUE)
  endif (CORE_INCLUDE_DIR AND CORE_LIBRARIES)

  find_path(CORE_INCLUDE_DIR NAMES CORE.h DOC "The directory containing the CORE include files")

  find_library(CORE_LIBRARIES NAMES core++ DOC "Path to the core++ library")
  
  get_filename_component(CORE_LIBRARIES_DIR ${CORE_LIBRARIES} PATH)

  FIND_PACKAGE_HANDLE_STANDARD_ARGS(CORE "DEFAULT_MSG" CORE_LIBRARIES CORE_INCLUDE_DIR )

endif(GMP_FOUND)
