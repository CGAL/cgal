#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Timer.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Epic;
typedef CGAL::Exact_predicates_exact_constructions_kernel       Epec;

template <typename K>
int
test_self_intersections(const std::string filename, const bool expected)
{
  typedef CGAL::Polyhedron_3<K>                                      Polyhedron;
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor  face_descriptor;

  std::ifstream input(filename);
  Polyhedron poly;

  if ( !input || !(input >> poly) ) {
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
    CGAL::parameters::vertex_index_map(get(CGAL::vertex_point, poly)));
  bool intersecting_1 = !intersected_tris.empty();

  std::cout << "Self-intersection test took " << timer.time() << " sec." << std::endl;
  std::cout << intersected_tris.size() << " pairs of triangles intersect." << std::endl;

  timer.reset();
  bool intersecting_2 = CGAL::Polygon_mesh_processing::does_self_intersect(poly,
    CGAL::parameters::vertex_index_map(get(CGAL::vertex_point, poly)));

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
  // If file(s) are provided, the associated expected result must also be provided.
  // Note that this expected value is a Boolean that is passed in command line
  // with either 'true' or 'false' (but not integers), that is for example:
  // > self_intersection_polyhedron_test data/U.off false

  // First test ----------------------------------------------------------------
  bool expected = false;
  std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");
  if(argc > 1) {
    assert(argc > 2);
    std::stringstream ss(argv[2]);
    ss >> std::boolalpha >> expected;
    assert(!ss.fail()); // make sure that argv[2] is either 'true' or 'false'
  }

  std::cout << "First test (Epic):" << std::endl;
  int r = test_self_intersections<Epic>(filename, expected);

  std::cout << "First test (Epec):" << std::endl;
  r += test_self_intersections<Epec>(filename, expected);

  // Second test ---------------------------------------------------------------
  expected = true;
  filename = (argc > 3) ? argv[3] : CGAL::data_file_path("meshes/mannequin-devil.off");
  if(argc > 3) {
    assert(argc > 4);
    std::stringstream ss(argv[4]);
    ss >> std::boolalpha >> expected;
    assert(!ss.fail());
  }

  std::cout << "Second test (Epic):" << std::endl;
  r += test_self_intersections<Epic>(filename, expected);

  std::cout << "Second test (Epec):" << std::endl;
  r += test_self_intersections<Epec>(filename, expected);

  // Third test ----------------------------------------------------------------
  expected = true;
  filename = (argc > 5) ? argv[5] : "data/overlapping_triangles.off";
  if(argc > 5) {
    assert(argc > 6);
    std::stringstream ss(argv[6]);
    ss >> std::boolalpha >> expected;
    assert(!ss.fail());
  }

  std::cout << "Third test (Epic):" << std::endl;
  r += test_self_intersections<Epic>(filename, expected);

  std::cout << "Third test (Epec):" << std::endl;
  r += test_self_intersections<Epec>(filename, expected);

  return r;
}
