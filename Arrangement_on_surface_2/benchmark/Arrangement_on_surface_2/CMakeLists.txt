cmake_minimum_required(VERSION 3.12...3.31)
project(Arrangement_on_surface_2_benchmarks)
find_package(Boost CONFIG ${BOOST_MIN_VERSION} REQUIRED COMPONENTS program_options)
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)

if(TARGET CGAL::CGAL_Basic_viewer)
  add_executable(do_intersect do_intersect.cpp)
  target_link_libraries(do_intersect PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 CGAL::CGAL_Basic_viewer Boost::program_options)
endif()

# Polylines with the Epeck kernel
add_executable(bench_arr_polylines bench_arr_polylines.cpp)
target_link_libraries(bench_arr_polylines PRIVATE CGAL::CGAL Boost::program_options)

# Polylines with the Simple_cartesian<Exact_rational> kernel
add_executable(bench_arr_polylines_rational bench_arr_polylines.cpp)
target_compile_definitions(bench_arr_polylines_rational PRIVATE CGAL_BENCH_KERNEL=2)
target_link_libraries(bench_arr_polylines_rational PRIVATE CGAL::CGAL Boost::program_options)
