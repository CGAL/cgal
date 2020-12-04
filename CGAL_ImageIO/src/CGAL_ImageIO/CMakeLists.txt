message("Configuring libCGAL_ImageIO")

collect_cgal_library(CGAL_ImageIO "")

include(CGAL_SetupCGAL_ImageIODependencies)

if(CGAL_HEADER_ONLY)
  set(keyword "INTERFACE")
endif()

cgal_setup_cgal_imageio_dependencies(CGAL_ImageIO ${keyword})

if(COMMAND add_config_flag)
  set(CGAL_HAS_IMAGEIO TRUE)
  add_config_flag(CGAL_HAS_IMAGEIO)
endif()

if(ZLIB_FOUND)
  get_dependency_version(ZLIB)
endif()

message("libCGAL_ImageIO is configured")
