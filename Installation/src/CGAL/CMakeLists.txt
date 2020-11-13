message("Configuring libCGAL")

collect_cgal_library(CGAL "")

include(CGAL_SetupCGALDependencies)

if(CGAL_HEADER_ONLY)
  set(keyword "INTERFACE")
endif()

cgal_setup_cgal_dependencies(CGAL ${keyword})

if(NOT CGAL_HEADER_ONLY)
  target_compile_definitions(CGAL INTERFACE CGAL_NOT_HEADER_ONLY=1)
endif()

if(NOT CGAL_DISABLE_GMP)
  get_dependency_version(GMP)
  get_dependency_version(MPFR)
endif()

message(
  STATUS
    "USING BOOST_VERSION = '${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}'"
)

message("libCGAL is configured")
