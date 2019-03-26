include(CGAL_SetupCGAL_CoreDependencies)

message("Configuring libCGAL_Core")
if(CGAL_Core_FOUND)

  collect_cgal_library(CGAL_Core "")

  if(CGAL_HEADER_ONLY)
    set(keyword INTERFACE)
  else()
    set(keyword PUBLIC)
  endif()

  CGAL_setup_CGAL_Core_dependencies(CGAL_Core ${keyword})

  message("libCGAL_Core is configured")
endif()

