#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Self_intersection_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;


int main(int, char** argv) {
  std::ifstream input(argv[1]);
  Polyhedron poly;

  if ( !input || !(input >> poly) ){
    std::cerr << "Error: can not read file.";
    return 1;
  }
  
  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  CGAL::self_intersections<K>(poly, back_inserter(intersected_tris));
  bool intersecting_1 = !intersected_tris.empty();
  CGAL_assertion(intersecting_1);

  std::cerr << "Self-intersection test took " << timer.time() << " sec." << std::endl;
  std::cerr << intersected_tris.size() << " pairs of triangles are intersecting." << std::endl;

  timer.reset();
  bool intersecting_2 = CGAL::do_self_intersect<K>(poly);
  CGAL_assertion(intersecting_1 == intersecting_2);

  std::cerr << "Is self-intersection test took " << timer.time() << " sec." << std::endl;
  std::cerr << (intersecting_2 ? "There is a self-intersection." : "There is no self-intersection.") << std::endl;

  return 0;
}
