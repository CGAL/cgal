include(CGAL_SetupCGAL_CoreDependencies)

message("Configuring libCGAL_Core")
if(CGAL_Core_FOUND)

  collect_cgal_library(CGAL_Core "")

  if(CGAL_HEADER_ONLY)
    set(keyword INTERFACE)
  else()
    set(keyword PUBLIC)
  endif()

  cgal_setup_cgal_core_dependencies(CGAL_Core ${keyword})

  message("libCGAL_Core is configured")
endif()
