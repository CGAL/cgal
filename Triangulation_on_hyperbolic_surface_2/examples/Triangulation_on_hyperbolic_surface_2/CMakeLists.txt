cmake_minimum_required(VERSION 3.12...3.31)

project( Triangulation_on_hyperbolic_surface_2_Examples )

# CGAL and its components
find_package( CGAL REQUIRED )


create_single_source_cgal_program("Triangulation_on_hyperbolic_surface_2.cpp")

if (GMP_FOUND AND MPFR_FOUND)
  create_single_source_cgal_program("Triangulation_on_hyperbolic_surface_2.cpp")
else()
  message("NOTICE: Some examples require the GMP and MPFR libraries, and will not be compiled.")
endif()
