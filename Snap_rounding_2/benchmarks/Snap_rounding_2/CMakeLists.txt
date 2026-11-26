cmake_minimum_required(VERSION 3.12...4.4)
project(Snap_rounding_2_Benchmarks)

find_package(CGAL REQUIRED)

find_package(Boost COMPONENTS program_options)
if(Boost_PROGRAM_OPTIONS_FOUND)
  create_single_source_cgal_program("bench_hot_pixel_snap_rounding_2.cpp")
  target_link_libraries(bench_hot_pixel_snap_rounding_2 PRIVATE Boost::program_options)
else()
  message("NOTICE: The benchmarks require Boost Program Options, and will not be compiled.")
endif()
