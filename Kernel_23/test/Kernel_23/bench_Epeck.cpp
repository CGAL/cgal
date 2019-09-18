#define CGAL_NO_CDT_2_WARNING 1
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Random.h>

#include <benchmark/benchmark.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;

typedef CGAL::Exact_intersections_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point Point;

static void bench_CDT_2_with_Epeck(benchmark::State& state) {
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

BENCHMARK(bench_CDT_2_with_Epeck)->RangeMultiplier(2)->Ranges({{16, 16<<4}, {0, 1}})->ArgNames({"", "degenerate"})->Unit(benchmark::kMillisecond);
BENCHMARK_MAIN();
