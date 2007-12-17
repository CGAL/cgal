# Try to find the CORE library shipped with CGAL
# CGAL_CORE_FOUND - system has CORE lib
# CGAL_CORE_INCLUDE_DIR - the CORE include directory
# CGAL_CORE_LIBRARIES   - the CGALCORE libraries

# TODO: support Windows and MacOSX

# CORE needs GMP

include(FindPackageHandleStandardArgs)

find_package(GMP QUIET)
if(GMP_FOUND)

  if (CGAL_CORE_INCLUDE_DIR AND CGAL_CORE_LIBRARIES)
    set(CGAL_CORE_FIND_QUIETLY TRUE)
  endif()

  # Find CORE include folder
  find_path(CGAL_CORE_INCLUDE_DIR NAMES CORE.h 
            PATHS ${CGAL_SOURCE_DIR}/include/CGAL/CORE
            DOC "The directory containing the CORE include files shipped with CGAL"
           )

  if ( AUTO_LINK_ENABLED )
  
    set(CGAL_CORE_LIBRARIES "" )
    set(CGAL_CORE_LIBRARIES_DIR "${CGAL_BINARY_DIR}/lib" )
    
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(CGAL_CORE "DEFAULT_MSG" CGAL_CORE_INCLUDE_DIR )
    
  else()
  
    # We cannot search for the core++ library because it is not yet compiled
    # => hard code the name
    if (WIN32)
      set(CGAL_CORE_LIBRARIES ${CGAL_BINARY_DIR}/lib/core++.lib)
    else()
      if(BUILD_SHARED_LIBS)
        set(CGAL_CORE_LIBRARIES ${CGAL_BINARY_DIR}/lib/libcore++.so)
      else(BUILD_SHARED_LIBS)
        set(CGAL_CORE_LIBRARIES ${CGAL_BINARY_DIR}/lib/libcore++.a)
      endif(BUILD_SHARED_LIBS)
    endif()
    
    set(CGAL_CORE_LIBRARIES_DIR "${CGAL_BINARY_DIR}/lib" )
    
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(CGAL_CORE "DEFAULT_MSG" CGAL_CORE_INCLUDE_DIR CGAL_CORE_LIBRARIES )
    
  endif()

endif()
