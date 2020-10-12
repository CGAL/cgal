include(CGAL_SetupCGAL_CoreDependencies)

message("Configuring libCGAL_Core")
if(CGAL_Core_FOUND)

  collect_cgal_library(CGAL_Core "")

  CGAL_setup_CGAL_Core_dependencies(CGAL_Core)

  message("libCGAL_Core is configured")
endif()

