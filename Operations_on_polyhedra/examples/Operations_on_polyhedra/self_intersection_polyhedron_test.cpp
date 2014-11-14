#undef NDEBUG
#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Self_intersection_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Triangle_3 Triangle;
 typedef CGAL::Polyhedron_3<K> Polyhedron;
//typedef CGAL::Surface_mesh<K::Point_3> Polyhedron;
typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;


int main(int, char** argv) {
  std::ifstream input(argv[1]);
  Polyhedron poly;
#if 0
  typedef boost::property_map<Polyhedron, CGAL::vertex_point_t>::type Ppmap;
  Ppmap point = get(CGAL::vertex_point_t(),poly);
  Point p = point[*vertices(poly).first];

#else

  if ( !input || !(input >> poly) ){
    std::cerr << "Error: can not read file.";
    return 1;
  }
  
  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  bool intersecting_1 = CGAL::self_intersect<K>(poly, back_inserter(intersected_tris)).first;
  assert(intersecting_1 == !intersected_tris.empty());

  std::cerr << "Self-intersection test took " << timer.time() << " sec." << std::endl;
  std::cerr << intersected_tris.size() << " pair of triangles are intersecting." << std::endl;

  timer.reset();
  bool intersecting_2 = CGAL::self_intersect<K>(poly);
  assert(intersecting_1 == intersecting_2);

  std::cerr << "Is self-intersection test took " << timer.time() << " sec." << std::endl;
  std::cerr << (intersecting_2 ? "There is a self-intersection." : "There is no self-intersection.") << std::endl;
#endif
  return 0;
}
