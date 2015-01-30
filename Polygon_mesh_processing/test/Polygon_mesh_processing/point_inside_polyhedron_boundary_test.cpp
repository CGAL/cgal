
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_inside_polygon_mesh_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include "point_inside_helpers.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;

int main(int, char** argv)
{
  std::ifstream input(argv[1]);
  Polyhedron poly;

  if ( !input || !(input >> poly) || poly.empty() ){
    std::cerr << "Error: can not read file.";
    return 1;
  }

  std::vector<Point> points;
  std::vector<bool> on_boundary;
  generate_near_boundary(poly, points, on_boundary);
  test(poly, points, on_boundary);

  points.clear();
  const int nb_query = (int)1.e6;
  points.reserve(nb_query);
  random_points<Point>(poly, nb_query, back_inserter(points));
  test(poly, points);

  //test compilation of constructor from AABB_tree  
  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> FGTP;
  typedef CGAL::AABB_traits<K, FGTP>    AABB_traits;
  typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;

  AABB_tree tree(faces(poly).first, faces(poly).second, poly);
  CGAL::Point_inside_polygon_mesh<Polyhedron, K> inside_test(tree);

}
