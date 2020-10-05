message("Configuring libCGAL")

collect_cgal_library(CGAL "")

include(CGAL_SetupCGALDependencies)

set(keyword "INTERFACE")

CGAL_setup_CGAL_dependencies(CGAL ${keyword})

if(NOT CGAL_DISABLE_GMP)
  get_dependency_version(GMP)
  get_dependency_version(MPFR)
endif()

message( STATUS "USING BOOST_VERSION = '${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}'" )

message("libCGAL is configured")
