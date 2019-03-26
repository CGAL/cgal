include(CGAL_SetupCGAL_CoreDependencies)

message("Configuring libCGAL_Core")
if(CGAL_Core_FOUND)
  if(NOT TARGET Exact_rational_support)
    set(CGAL_Core_EXTRA_TARGET Exact_rational_support)
  endif()
  collect_cgal_library(CGAL_Core "")

  if(CGAL_HEADER_ONLY)
    set(keyword INTERFACE)
  else()
    set(keyword PUBLIC)
  endif()

  CGAL_setup_CGAL_Core_dependencies(CGAL_Core ${keyword})

  message("libCGAL_Core is configured")
endif()

