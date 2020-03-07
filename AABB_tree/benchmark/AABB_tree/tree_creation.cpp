#include <benchmark/benchmark.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef K::Segment_3 Segment;
typedef K::Point_3 Point_3;


Surface_mesh mesh;

static void BM_TreeCreation(benchmark::State& state)
{
  for (auto _ : state)
  {
    benchmark::DoNotOptimize([]() {
      Tree tree{mesh.faces_begin(), mesh.faces_end(), mesh};
      tree.build();
      return 0;
    }());
  }
}
BENCHMARK(BM_TreeCreation);

static void BM_Intersections(benchmark::State& state)
{
  Point_3 p(-0.5, 0.03, 0.04);
  Point_3 q(-0.5, 0.04, 0.06);

  Tree tree{mesh.faces_begin(), mesh.faces_end(), mesh};
  tree.accelerate_distance_queries();

  Segment segment_query(p, q);
  for (auto _ : state)
  {
    benchmark::DoNotOptimize([&]() {
      tree.number_of_intersected_primitives(segment_query);
      Point_3 point_query(2.0, 2.0, 2.0);
      Point_3 closest = tree.closest_point(point_query);
      return 0;
    }());
  }
}
BENCHMARK(BM_Intersections);


int main(int argc, char** argv)
{
  const char* default_file = "data/handle.off";
  const char* filename = argc > 2? argv[2] : default_file;

  {
    std::ifstream input(filename);
    input >> mesh;
  }

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();

  return EXIT_SUCCESS;
}
