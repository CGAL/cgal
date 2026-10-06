cmake_minimum_required(VERSION 3.12...3.31)
project(do_intersect)
find_package(Boost CONFIG ${BOOST_MIN_VERSION} REQUIRED COMPONENTS program_options)
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt6)
add_executable(do_intersect do_intersect.cpp)
target_link_libraries(do_intersect PRIVATE CGAL::CGAL CGAL::CGAL_Qt6 CGAL::CGAL_Basic_viewer Boost::program_options)

# Performance of the operations that go through the arrangement geometry
# traits; used, e.g., to compare builds before and after a change.
add_executable(arr_traits_benchmark arr_traits_benchmark.cpp)
target_link_libraries(arr_traits_benchmark PRIVATE CGAL::CGAL)
