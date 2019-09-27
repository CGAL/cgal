#define CGAL_NO_CDT_2_WARNING 1
// #define CGAL_NO_STATIC_FILTERS 1
// #define CGAL_PROFILE 1
#include <CGAL/intersections.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Random.h>

// surface mesh
#include <CGAL/Polyhedron_3.h>

// nef
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>

#include <benchmark/benchmark.h>

#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>

template <typename K>
static void bench_nef(benchmark::State& state) {
  typedef CGAL::Polyhedron_3<K> Exact_polyhedron;

  typedef CGAL::Nef_polyhedron_3<K,
                                 CGAL::SNC_indexed_items,
                                 bool> Nef_polyhedron; 
  std::ifstream off_a("data/couplingdown.off");
  Exact_polyhedron poly;
  off_a >> poly;
  if(!off_a) state.SkipWithError("Failed to read \"data/couplingdown.off\"!");
  const Nef_polyhedron nef_a{poly};
  off_a.close();
  poly.clear();
  std::ifstream off_b("data/elephant.off");
  off_b >> poly;
  if(!off_a) state.SkipWithError("Failed to read \"data/elephant.off\"!");
  const Nef_polyhedron nef_b{poly};
  off_b.close();
  poly.clear();
  for(auto _ : state) {
    Nef_polyhedron result = nef_b - nef_a;
    benchmark::DoNotOptimize(result);
  }
}
template <typename K>
static void bench_CDT_2(benchmark::State& state) {
  typedef CGAL::Exact_intersections_tag Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
  typedef typename CDT::Point Point;
  CGAL::Random random{0};
  for(auto _ : state) {
    CDT cdt;
    const int nb_of_segments_x = state.range(0);
    const int nb_of_segments_y = state.range(0);
    const int non_degenerate = 1-state.range(1);
    const std::size_t nb_of_vertices = nb_of_segments_x * nb_of_segments_y +
      2 * (nb_of_segments_x + nb_of_segments_y);
    for(int i = 0; i < nb_of_segments_x; ++i) {
      cdt.insert_constraint(Point(i+non_degenerate*random.get_double(0., 0.5),
                                  -1),
                            Point(i+non_degenerate*random.get_double(0., 0.5),
                                  nb_of_segments_y + 1.));
    }
    for(int i = 0; i < nb_of_segments_y; ++i) {
      cdt.insert_constraint(Point(-1,
                                  i+non_degenerate*random.get_double(0., 0.5)),
                            Point(nb_of_segments_x + 1.,
                                  i+non_degenerate*random.get_double(0., 0.5)));
    }
    if(cdt.number_of_vertices() != nb_of_vertices) {
      state.SkipWithError("Wrong CDT");
    }
  }
}

BENCHMARK_TEMPLATE(bench_nef, CGAL::Epeck)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(bench_nef, CGAL::Atomic_ref_counted_epeck)->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(bench_nef, CGAL::Thread_safe_epeck)->Unit(benchmark::kMillisecond);

BENCHMARK_TEMPLATE(bench_CDT_2, CGAL::Epeck)->RangeMultiplier(2)->Ranges({{16, 16<<4}, {0, 1}})->ArgNames({"", "degenerate"})->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(bench_CDT_2, CGAL::Atomic_ref_counted_epeck)->RangeMultiplier(2)->Ranges({{16, 16<<4}, {0, 1}})->ArgNames({"", "degenerate"})->Unit(benchmark::kMillisecond);
BENCHMARK_TEMPLATE(bench_CDT_2, CGAL::Thread_safe_epeck)->RangeMultiplier(2)->Ranges({{16, 16<<4}, {0, 1}})->ArgNames({"", "degenerate"})->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
