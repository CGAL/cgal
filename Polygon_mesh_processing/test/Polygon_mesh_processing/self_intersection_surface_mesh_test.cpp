#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel       Epec;

template <typename K>
int
test_self_intersections(const char* filename, const bool expected)
{
  typedef CGAL::Surface_mesh<typename K::Point_3>                Mesh;
  typedef typename boost::graph_traits<Mesh>::face_descriptor    face_descriptor;

  std::ifstream input(filename);
  Mesh m;

  if ( !input || !(input >> m) ){
    std::cerr << "Error: can not read file.";
    return 1;
  }

  std::cout << "Reading file: " << filename << std::endl;

  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  CGAL::Polygon_mesh_processing::self_intersections(
    m,
    std::back_inserter(intersected_tris),
    CGAL::Polygon_mesh_processing::parameters::vertex_index_map(get(CGAL::vertex_point, m)));
  bool intersecting_1 = !intersected_tris.empty();

  std::cout << "self_intersections test took " << timer.time() << " sec." << std::endl;
  std::cout << intersected_tris.size() << " pairs of triangles are intersecting." << std::endl;

  timer.reset();
  bool intersecting_2 = CGAL::Polygon_mesh_processing::does_self_intersect(m,
    CGAL::Polygon_mesh_processing::parameters::vertex_index_map(get(CGAL::vertex_point, m)));

  std::cerr << "does_self_intersect test took " << timer.time() << " sec." << std::endl;
  std::cerr << (intersecting_2 ? "There is a self-intersection." :
                                 "There are no self-intersections.") << std::endl;

  assert(intersecting_1 == intersecting_2);
  assert(intersecting_1 == expected);

  std::cout << filename << " passed the tests." << std::endl << std::endl;

  return 0;
}

int main(int argc, char** argv)
{
  const char* filename_ele = (argc > 1) ? argv[1] : "data/elephant.off";
  const char* filename_man = (argc > 2) ? argv[2] : "data/mannequin-devil.off";
  const char* filename_ove = (argc > 3) ? argv[3] : "data/overlapping_triangles.off";

  std::cout << "Testing with Epic:" << std::endl;
  int r = test_self_intersections<Epic>(filename_ele, false);
  r += test_self_intersections<Epic>(filename_man, true);
  r += test_self_intersections<Epic>(filename_ove, true);

  std::cout << "Testing with Epec:" << std::endl;
  r += test_self_intersections<Epec>(filename_ele, false);
  r += test_self_intersections<Epec>(filename_man, true);
  r += test_self_intersections<Epec>(filename_ove, true);

  return r;
}
