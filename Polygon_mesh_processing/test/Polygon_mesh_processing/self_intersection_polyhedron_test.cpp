#include <cstdlib>
#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel       Epec;

template <typename K>
int
test_self_intersections(const char* filename, const bool expected)
{
  typedef CGAL::Polyhedron_3<K>                                      Polyhedron;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor  face_descriptor;

  std::ifstream input(filename);
  Polyhedron poly;

  if ( !input || !(input >> poly) ){
    std::cerr << "Error: cannot read file: " << filename << std::endl;
    return 1;
  }

  std::cout << "Reading file: " << filename << std::endl;

  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  CGAL::Polygon_mesh_processing::self_intersections(
    poly,
    std::back_inserter(intersected_tris),
    CGAL::Polygon_mesh_processing::parameters::vertex_index_map(get(CGAL::vertex_point, poly)));
  bool intersecting_1 = !intersected_tris.empty();

  std::cout << "Self-intersection test took " << timer.time() << " sec." << std::endl;
  std::cout << intersected_tris.size() << " pairs of triangles intersect." << std::endl;

  timer.reset();
  bool intersecting_2 = CGAL::Polygon_mesh_processing::does_self_intersect(poly,
      CGAL::Polygon_mesh_processing::parameters::vertex_index_map(get(CGAL::vertex_point, poly)));

  std::cout << "does_self_intersect test took " << timer.time() << " sec." << std::endl;
  std::cout << (intersecting_2 ? "There is a self-intersection." :
                                 "There are no self-intersections.") << std::endl;

  assert(intersecting_1 == intersecting_2);
  assert(intersecting_1 == expected);

  std::cout << filename << " passed the tests." << std::endl << std::endl;

  return 0;
}

int main(int argc, char** argv)
{
  const char* first_file_name = (argc > 1) ? argv[1] : "data/elephant.off";
  if(argc > 1) assert(argc > 2);
  const bool first_expected = (argc > 1) ? std::atoi(argv[2]) : false;

  const char* second_file_name = (argc > 3) ? argv[3] : "data/mannequin-devil.off";
  if(argc > 3) assert(argc > 4);
  const bool second_expected = (argc > 3) ? std::atoi(argv[4]) : true;

  const char* third_file_name = (argc > 5) ? argv[5] : "data/overlapping_triangles.off";
  if(argc > 5) assert(argc > 6);
  const bool third_expected = (argc > 5) ? std::atoi(argv[6]) : true;


  std::cout << "Testing with Epic:" << std::endl;
  int r = test_self_intersections<Epic>(first_file_name, first_expected);
  r += test_self_intersections<Epic>(second_file_name, second_expected);
  r += test_self_intersections<Epic>(third_file_name, third_expected);

  std::cout << "Testing with Epec:" << std::endl;
  r += test_self_intersections<Epec>(first_file_name, first_expected);
  r += test_self_intersections<Epec>(second_file_name, second_expected);
  r += test_self_intersections<Epec>(third_file_name, third_expected);

  return r;
}
