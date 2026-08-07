# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.12...3.31)
project(Triangulation_on_hyperbolic_surface_2_Tests)

find_package(CGAL REQUIRED)

create_single_source_cgal_program("hs_test_circular_kernel.cpp")
create_single_source_cgal_program("hs_test_complex.cpp")
create_single_source_cgal_program("hs_test_domain.cpp")
create_single_source_cgal_program("hs_test_factory.cpp")
create_single_source_cgal_program("hs_test_isometry.cpp")
create_single_source_cgal_program("hs_test_lazy_exact_nt.cpp")
create_single_source_cgal_program("hs_test_triangulation.cpp")

if (GMP_FOUND AND MPFR_FOUND)
  create_single_source_cgal_program("hs_test_Delaunay_triangulation.cpp")
  create_single_source_cgal_program("hs_test_espilon_net.cpp")
else()
  message("NOTICE: Some tests require the GMP and MPFR libraries, and will not be compiled.")
endif()
