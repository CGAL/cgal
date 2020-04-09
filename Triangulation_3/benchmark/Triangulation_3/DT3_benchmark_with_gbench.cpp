#include <CGAL/Real_timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Random.h>

#include <iostream>
#include <fstream>
#include <benchmark/benchmark.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::Point_3                                           Point_3;

typedef CGAL::Triangulation_data_structure_3<
                  CGAL::Triangulation_vertex_base_3<K>,
              CGAL::Triangulation_cell_base_3<K>,
                  CGAL::Parallel_tag>                                Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>               PDT;

// global variables used by bench_dt3
int argc;
char** argv;



void bench_dt3(benchmark::State& state) {
  CGAL::get_default_random() = CGAL::Random(0);

  std::vector<Point_3> points;
  Point_3 p;

  std::ifstream in(argv[1]);
  while(in >> p)
    points.push_back(p);

  for(auto _ : state) {
    CGAL::Bbox_3 bb = CGAL::bounding_box(points.begin(), points.end()).bbox();
    PDT::Lock_data_structure locking_ds(bb, 50);

    PDT pdt(points.begin(), points.end(), &locking_ds);
  }
  return;
}
BENCHMARK(bench_dt3)->Unit(benchmark::kMillisecond);;


int main(int argc, char* argv[])
{
  benchmark::Initialize(&argc, argv);
  ::argc = argc;
  ::argv = argv;
  benchmark::RunSpecifiedBenchmarks();
}
