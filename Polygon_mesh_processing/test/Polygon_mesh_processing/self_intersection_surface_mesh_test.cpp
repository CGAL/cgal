#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;


int main(int, char** argv) {
  std::ifstream input(argv[1]);
  Mesh m;

  if ( !input || !(input >> m) ){
    std::cerr << "Error: can not read file.";
    return 1;
  }
  
  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  CGAL::Polygon_mesh_processing::self_intersections<K>
    (m, back_inserter(intersected_tris));
  bool intersecting_1 = !intersected_tris.empty();

  std::cerr << "self_intersections test took " << timer.time() << " sec." << std::endl;
  std::cerr << intersected_tris.size() << " pairs of triangles are intersecting." << std::endl;

  timer.reset();
  bool intersecting_2
    = CGAL::Polygon_mesh_processing::is_self_intersecting<K>(m);
  CGAL_assertion(intersecting_1 == intersecting_2);

  std::cerr << "is_self_intersecting test took " << timer.time() << " sec." << std::endl;
  std::cerr << (intersecting_2 ? "There is a self-intersection." : "There is no self-intersection.") << std::endl;

  return 0;
}
