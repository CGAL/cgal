
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include "point_inside_helpers.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;

int main(int argc, char** argv)
{
  const char* filename = (argc > 1) ? argv[1] : "data/elephant.off";
  std::ifstream input(filename);
  Polyhedron poly;

  if ( !input || !(input >> poly) || poly.empty() ){
    std::cerr << "Error: can not read file.";
    return 1;
  }

  std::vector<Point> points;
  std::vector<bool> on_boundary;
  generate_near_boundary(poly, points, on_boundary);
  inside_test(poly, points, on_boundary);

  points.clear();
  const int nb_query = (int)1.e4;
  points.reserve(nb_query);
  random_points<Point>(poly, nb_query, back_inserter(points));
  inside_test(poly, points);

  //test compilation of constructor from AABB_tree
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> FGTP;
  typedef CGAL::AABB_traits<K, FGTP>    AABB_traits;
  typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;

  AABB_tree tree(faces(poly).first, faces(poly).second, poly);
  CGAL::Side_of_triangle_mesh<Polyhedron, K> inside_tester(tree);

  CGAL::Bounded_side bs = inside_tester(CGAL::ORIGIN);
  std::cout << "Origin is " << bs << std::endl;

  return 0;
}
