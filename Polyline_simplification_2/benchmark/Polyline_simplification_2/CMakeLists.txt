cmake_minimum_required(VERSION 3.12...3.29)
project(Polyline_simplification_2_Benchmarks)

# Add Google Benchmark
find_package(benchmark REQUIRED)

# Add CGAL
find_package(CGAL REQUIRED)

add_executable(benchmark_simplify benchmark_simplify.cpp)
target_link_libraries(benchmark_simplify PRIVATE benchmark::benchmark CGAL::CGAL CGAL::Data)
