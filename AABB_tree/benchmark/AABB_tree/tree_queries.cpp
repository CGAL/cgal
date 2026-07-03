#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_mesh_io.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>

#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>

namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Simple_cartesian<double>;

using Point   = K::Point_3;
using Vector  = K::Vector_3;
using Segment = K::Segment_3;
using Ray     = K::Ray_3;
using Line    = K::Line_3;
using Plane   = K::Plane_3;

using Mesh = CGAL::Surface_mesh<Point>;

using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits    = CGAL::AABB_traits_3<K, Primitive>;
using Tree      = CGAL::AABB_tree<Traits>;

CGAL::Random rng;

template<class Kernel = K>
typename Kernel::Point_3 random_point(const CGAL::Bbox_3& bb){
  return typename Kernel::Point_3(rng.get_double(bb.xmin(), bb.xmax()),
                                  rng.get_double(bb.ymin(), bb.ymax()),
                                  rng.get_double(bb.zmin(), bb.zmax()));
}

Vector random_vector(){
  Vector v;
  do{
    v = Vector(rng.get_double(-1., 1.),
               rng.get_double(-1., 1.),
               rng.get_double(-1., 1.));
  } while(v.squared_length() < 1e-12);
  return v;
}

template<class QueryGenerator, class Function>
double benchmark(QueryGenerator gen,
                 Function f,
                 std::size_t nb_queries = 500000)
{
  CGAL::Real_timer timer;
  timer.start();
  for(std::size_t i=0; i<nb_queries; ++i)
    f(gen());
  timer.stop();
  return double(nb_queries) / timer.time();
}

void benchmark_queries(const Mesh& mesh)
{
  Tree tree(faces(mesh).first, faces(mesh).second, mesh);
  tree.build();

  CGAL::Bbox_3 bb = PMP::bbox(mesh);
  constexpr std::size_t N = 100000;
  // Warm-up
  for(int i=0;i<10000;++i)
    tree.do_intersect(Segment(random_point(bb), random_point(bb)));
  std::cout << std::fixed << std::setprecision(0);

  auto print = [](const std::string& name, double value)  {
    std::cout << std::setw(40) << std::left
              << name
              << value
              << " queries/s\n";
  };

#define BENCHMARK(QueryType, Generator)                                      \
  {                                                                          \
    std::cout << "\n=== " #QueryType " ===\n";                               \
                                                                             \
    print("do_intersect",                                                    \
      benchmark(Generator,                                                   \
      [&](const QueryType& q)                                                \
      { tree.do_intersect(q); }, N));                                        \
                                                                             \
    print("any_intersected_primitive",                                       \
      benchmark(Generator,                                                   \
      [&](const QueryType& q)                                                \
      { tree.any_intersected_primitive(q); }, N));                           \
                                                                             \
    print("any_intersection",                                                \
      benchmark(Generator,                                                   \
      [&](const QueryType& q)                                                \
      { tree.any_intersection(q); }, N));                                    \
                                                                             \
    print("number_of_intersected_primitives",                                \
      benchmark(Generator,                                                   \
      [&](const QueryType& q)                                                \
      { tree.number_of_intersected_primitives(q); }, N));                    \
                                                                             \
    print("all_intersected_primitives",                                      \
      benchmark(Generator,                                                   \
      [&](const QueryType& q)                                                \
      {                                                                      \
        std::vector<Tree::Primitive_id> ids;                                 \
        tree.all_intersected_primitives(q, std::back_inserter(ids));          \
      }, N));                                                                \
                                                                             \
    print("all_intersections",                                               \
      benchmark(Generator,                                                   \
      [&](const QueryType& q)                                                \
      {                                                                      \
        using Intersection =                                                 \
          typename Tree::template Intersection_and_primitive_id<QueryType>::Type; \
                                                                             \
        std::vector<Intersection> out;                                       \
        tree.all_intersections(q, std::back_inserter(out));                  \
      }, N));                                                                \
  }

  BENCHMARK(Segment, [&](){ return Segment(random_point(bb), random_point(bb)); });
  BENCHMARK(Ray, [&](){ return Ray(random_point(bb), random_vector()); });
  BENCHMARK(Line, [&](){ return Line(random_point(bb), random_point(bb)); });
  BENCHMARK(Plane, [&](){ return Plane(random_point(bb), random_vector()); });
#undef BENCHMARK
}

template<class Kernel>
void benchmark_kernel(const std::string &filename, const std::string &Kernel_name)
{
  using Segment = typename Kernel::Segment_3;

  using Mesh = CGAL::Surface_mesh<typename Kernel::Point_3>;
  using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
  using Traits    = CGAL::AABB_traits_3<Kernel, Primitive>;
  using Tree      = CGAL::AABB_tree<Traits>;
  const int N = 100000;

  Mesh mesh;
  CGAL::IO::read_polygon_mesh(filename, mesh);
  Tree tree(faces(mesh).first, faces(mesh).second, mesh);
  tree.build();
  CGAL::Bbox_3 bb = PMP::bbox(mesh);

  CGAL::Real_timer t;
  t.start();
  for(std::size_t i=0;i<N;++i){
      Segment s(random_point<Kernel>(bb), random_point<Kernel>(bb));
      using Intersection = typename Tree::template Intersection_and_primitive_id<Segment>::Type;
      std::vector<Intersection> out;
      tree.all_intersections(s, std::back_inserter(out));
  }
  t.stop();

  std::cout << std::setw(40)
            << Kernel_name
            << N/t.time()
            << '\n';
}

void benchmark_distances(const Mesh& mesh, std::size_t N = 300000)
{
  Tree tree(faces(mesh).first, faces(mesh).second, mesh);
  tree.build();
  tree.accelerate_distance_queries();

  const CGAL::Bbox_3 bb = PMP::bbox(mesh);

  std::vector<Point> queries;
  queries.reserve(N);
  for(std::size_t i=0; i<N; ++i)
    queries.push_back(random_point(bb));

  auto benchmark = [&](const std::string& name, auto&& f)  {
    CGAL::Real_timer timer;
    timer.start();
    for(const Point& p : queries)
      f(p);
    timer.stop();

    std::cout << std::setw(32) << std::left
              << name
              << std::fixed << std::setprecision(0)
              << N / timer.time()
              << " queries/s\n";
  };

  benchmark("closest_point()", [&](const Point& p){ tree.closest_point(p); });
  benchmark("squared_distance()", [&](const Point& p){ tree.squared_distance(p); });
  benchmark("closest_point_and_primitive()", [&](const Point& p){ tree.closest_point_and_primitive(p); });
}

int main(int argc, char** argv)
{
  if(argc != 2){
    std::cerr << "Usage: " << argv[0] << " mesh.off\n";
    return EXIT_FAILURE;
  }

  Mesh mesh;
  std::string filename = argv[1];

  if(!CGAL::IO::read_polygon_mesh(filename, mesh))  {
    std::cerr << "Cannot read mesh\n";
    return EXIT_FAILURE;
  }
  // benchmark_queries(mesh);
  // benchmark_kernel<CGAL::Simple_cartesian<double>>(filename, "Simple_cartesian<double>");
  // benchmark_kernel<CGAL::Simple_cartesian<float>>(filename, "Simple_cartesian<float>");
  // benchmark_kernel<CGAL::Cartesian<double>>(filename, "Cartesian<double>");
  // benchmark_kernel<CGAL::Cartesian<float>>(filename, "Cartesian<float>>");
  // benchmark_kernel<CGAL::Epick>(filename, "Epick");
  benchmark_distances(mesh);

  return EXIT_SUCCESS;
}