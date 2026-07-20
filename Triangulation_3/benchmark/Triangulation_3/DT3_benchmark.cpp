//#define CGAL_PROFILE
//#define CGAL_USE_SSE2_FABS
//#define CGAL_USE_SSE2_MAX
//#define CGAL_MSVC_USE_STD_FABS  // use this one with precise

#define INDEX_STORAGE 1

#define CGAL_NDEBUG 1
#define NDEBUG 1

#if PARALLEL
#  include "CGAL/Bbox_3.h"
#endif
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "CGAL/Triangulation_data_structure_3.h"
#include "CGAL/Triangulation_vertex_base_3.h"


#include <CGAL/Delaunay_triangulation_3.h>

#include <iostream>
#include <string>
#include <fstream>
#include <locale>

#include <benchmark/benchmark.h>

#if PARALLEL
using Concurrent_tag = CGAL::Parallel_tag;
#else
using Concurrent_tag = CGAL::Sequential_tag;
#endif
typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
#ifdef  INDEX_STORAGE
typedef CGAL::VertexWithPoint<K> Vb;
typedef CGAL::Cell4Delaunay<K> Cb;
typedef CGAL::Index_tag Tds_type_tag;
#else
typedef CGAL::Triangulation_vertex_base_3<K> Vb;
typedef CGAL::Delaunay_triangulation_cell_base_3<K> Cb;
typedef CGAL::Handle_tag Tds_type_tag;
#endif
typedef CGAL::Triangulation_data_structure_3<Vb, Cb, Concurrent_tag, Tds_type_tag> Tds;
typedef CGAL::Delaunay_triangulation_3<K,Tds>                DT;
typedef DT::Point                                            Point_3;

// global variables used by bench_dt3
int argc;
char** argv;



void bench_dt3(benchmark::State& state) {
  std::locale loc = std::locale()
      .combine<std::numpunct<char>>(std::locale("en_US.UTF8"));
  std::cout.imbue(loc);

  int M = 100; // Number of times to compute the triangulation

  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/ocean_r.xyz");
  if(argc > 2) {
    auto M_ = std::atoi(argv[2]);
    if(M_ <= 0) {
      std::cerr << "Invalid number of iterations: " << M_ << ". Using default value of 100." << std::endl;
    } else {
      M = M_;
    }
  }
  std::ifstream in(filename.c_str());
  std::vector<Point_3> points;
  Point_3 p, q;

  while(in >> p ){
    points.push_back(p);
  }

  for(auto _ : state) {
  std::cout << "Compute triangulation "<< M << " times" << std::endl;
  for(int i = 0; i < M; i++){

    DT dt;
    dt.insert(points.begin(), points.end());
    std::cout << dt.number_of_cells() << std::endl;
  }
  }
}

BENCHMARK(bench_dt3)->Unit(benchmark::kMillisecond);


int main(int argc, char* argv[])
{
  benchmark::Initialize(&argc, argv);
  ::argc = argc;
  ::argv = argv;
  benchmark::RunSpecifiedBenchmarks();
}
