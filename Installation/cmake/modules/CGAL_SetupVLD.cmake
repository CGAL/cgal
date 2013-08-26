function(CGAL_setupVLD)

  if(DEFINED ENV{CGAL_ENABLE_VLD} AND NOT DEFINED CGAL_ENABLE_VLD)
    set(CGAL_ENABLE_VLD $ENV{CGAL_ENABLE_VLD})
  endif()

  if(CGAL_ENABLE_VLD)
    find_path(VLD_INCLUDE_DIR vld.h
      HINTS $ENV{VLD_INCLUDE_DIR}
            $ENV{VLD_HOME}/include
            $ENV{VLD_HOME}
      DOC "Path of the Visual Leak Detector header vld.h"
    )

    find_library(VLD_LIBRARY
      NAMES vld-x86 vld
      HINTS ${VLD_LIBRARY_DIR}
            $ENV{VLD_LIBRARY_DIR}
            $ENV{VLD_HOME}/lib
            $ENV{VLD_HOME}
            ${VLD_INCLUDE_DIR}
            ${VLD_INCLUDE_DIR}/lib
      DOC "Path of the Visual Leak Detector library vld.lib"
    )
    if(VLD_LIBRARY)
      get_filename_component(VLD_LIBRARY_DIR "${VLD_LIBRARY}" PATH CACHE)
    endif()

    if(VLD_INCLUDE_DIR AND VLD_LIBRARY_DIR)
      set(VLD_FOUND 1 PARENT_SCOPE)
      include_directories(${VLD_INCLUDE_DIR})
      link_directories(${VLD_LIBRARY_DIR})
      add_definitions(-DCGAL_ENABLE_VLD)
    endif()

  endif(CGAL_ENABLE_VLD)
endfunction()
