# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


cmake_minimum_required(VERSION 3.1...3.15)
project( Arrangement_on_surface_2_Tests )

enable_testing()


find_package(CGAL QUIET COMPONENTS Core)

if ( CGAL_FOUND )

  include(${CGAL_USE_FILE})

  if(COMMAND target_compile_options)
    # Since CMake-2.8.12: New CMake script, that defines the targets and
    # the CTest test cases.
    include( ${CMAKE_CURRENT_SOURCE_DIR}/cgal_test.cmake )
  else()
    # If CMake version is <= 2.8.11, use the usual CMake script.

    # create a target per cppfile
    file(GLOB cppfiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

    create_single_source_cgal_program("construction_test_suite_generator.cpp")
    create_single_source_cgal_program("ex_kernel_point.cpp")
    create_single_source_cgal_program("ex_kernel_segment.cpp")
    create_single_source_cgal_program("test_arc_polycurve.cpp")
    create_single_source_cgal_program("test_batched_point_location.cpp")
    create_single_source_cgal_program("test_bezier_polycurve.cpp")
    create_single_source_cgal_program("test_circular_arc_polycurve.cpp")
    create_single_source_cgal_program("test_conic_polycurve.cpp")
    create_single_source_cgal_program("test_construction.cpp")
    create_single_source_cgal_program("test_data_traits.cpp")
    create_single_source_cgal_program("test_do_equal.cpp")
    create_single_source_cgal_program("test_do_intersect.cpp")
    create_single_source_cgal_program("test_dual.cpp")
    create_single_source_cgal_program("test_insertion.cpp")
    create_single_source_cgal_program("test_io.cpp")
    create_single_source_cgal_program("test_iso_verts.cpp")
    create_single_source_cgal_program("test_observer.cpp")
    create_single_source_cgal_program("test_overlay.cpp")
    create_single_source_cgal_program("test_point_location.")
    create_single_source_cgal_program("test_point_location_dynamic.cpp")
    create_single_source_cgal_program("test_polycurve_intersection.cpp")
    create_single_source_cgal_program("test_polylines.cpp")
    create_single_source_cgal_program("test_rational_function_traits_2.cpp")
    create_single_source_cgal_program("test_removal.cpp")
    create_single_source_cgal_program("test_sgm.cpp")
    create_single_source_cgal_program("test_spherical_merge.cpp")
    create_single_source_cgal_program("test_spherical_overlay.cpp")
    create_single_source_cgal_program("test_spherical_removal.cpp")
    create_single_source_cgal_program("test_tags.cpp")
    create_single_source_cgal_program("test_traits_adaptor.cpp")
    create_single_source_cgal_program("test_traits_dispatching.cpp")
    create_single_source_cgal_program("test_unbounded_rational_direct_insertion.cpp")
    create_single_source_cgal_program("test_unbounded_rational_insertion.cpp")
    create_single_source_cgal_program("test_unbounded_removal.cpp")
    create_single_source_cgal_program("test_vertical_decomposition.cpp")
    create_single_source_cgal_program("test_vert_ray_shoot_vert_segments.cpp")
    create_single_source_cgal_program("test_zone.cpp")
    if(NOT CGAL_DISABLE_GMP)
      create_single_source_cgal_program("test_traits.cpp")
    else()
      message(STATUS "NOTICE: Some tests require the CGAL_Core library, and will not be compiled.")
    endif()
  endif()

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
