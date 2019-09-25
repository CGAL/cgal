#define CGAL_NO_CDT_2_WARNING 1
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

typedef CGAL::Exact_predicates_exact_constructions_kernel K;

typedef CGAL::Exact_intersections_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, Itag> CDT;
typedef CDT::Point Point;

typedef CGAL::Polyhedron_3<K> Exact_polyhedron;

typedef CGAL::Nef_polyhedron_3<K,
			       CGAL::SNC_indexed_items,
			       bool> Nef_polyhedron;

Nef_polyhedron nef_a, nef_b;

class MyNefFixture : public benchmark::Fixture {
public:
  void SetUp(const ::benchmark::State&) {
    std::ifstream off_a("data/couplingdown.off");
    Exact_polyhedron poly;
    off_a >> poly;
    nef_a = Nef_polyhedron(poly);
    off_a.close();
    poly.clear();
    std::ifstream off_b("data/elephant.off");
    off_b >> poly;
    nef_b = Nef_polyhedron(poly);
  }

  void TearDown(const ::benchmark::State&) {
    nef_a.clear();
    nef_b.clear();
  }
};

BENCHMARK_DEFINE_F(MyNefFixture, bench_nef_with_Epeck)(benchmark::State& state) {
  for(auto _ : state) {
    Nef_polyhedron result = nef_b - nef_a;
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK_REGISTER_F(MyNefFixture, bench_nef_with_Epeck)->Unit(benchmark::kMillisecond);

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
