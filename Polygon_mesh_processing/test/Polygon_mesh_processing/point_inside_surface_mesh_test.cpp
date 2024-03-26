#include <CGAL/Side_of_triangle_mesh.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/boost/graph/helpers.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include "point_inside_helpers.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel Epec;


template <class K>
int test_surface_mesh(const std::string filename)
{
  typedef typename K::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> Mesh;
  Mesh mesh;

  std::ifstream input(filename);
  if (!input || !(input >> mesh)){
    std::cerr << "Error: can not read file.";
    return 1;
  }

  std::vector<Point> points;
  std::vector<bool> on_boundary;
  generate_near_boundary(mesh, points, on_boundary);
  inside_test(mesh, points, on_boundary);

  points.clear();
  const int nb_query = (int)1.e4;
  points.reserve(nb_query);
  random_points<Point>(mesh, nb_query, back_inserter(points));
  inside_test(mesh, points);

  //test compilation of constructor from AABB_tree
  typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> FGTP;
  typedef CGAL::AABB_traits<K, FGTP>    AABB_traits;
  typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;

  AABB_tree tree(faces(mesh).first, faces(mesh).second, mesh);
  CGAL::Side_of_triangle_mesh<Mesh, K> inside_tester(tree);

  CGAL::Bounded_side bs = inside_tester(CGAL::ORIGIN);
  std::cout << "Origin is " << bs << std::endl;
  return 0;
}

int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  assert(test_surface_mesh<Epic>(filename) == 0);
  assert(test_surface_mesh<Epec>(filename) == 0);
  return 0;
}
