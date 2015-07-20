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


int main(int argc, char** argv)
{
  const char* filename = (argc > 1) ? argv[1] : "data/elephant.off";
  std::ifstream input(filename);
  Mesh m;

  if ( !input || !(input >> m) ){
    std::cerr << "Error: can not read file.";
    return 1;
  }
  
  CGAL::Timer timer;
  timer.start();

  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  CGAL::Polygon_mesh_processing::self_intersections(m,
    std::back_inserter(intersected_tris),
    CGAL::Polygon_mesh_processing::parameters::vertex_index_map(get(CGAL::vertex_point, m)));
  bool intersecting_1 = !intersected_tris.empty();
  assert(!intersecting_1);

  std::cerr << "self_intersections test took " << timer.time() << " sec." << std::endl;
  std::cerr << intersected_tris.size() << " pairs of triangles are intersecting." << std::endl;

  timer.reset();
  bool intersecting_2 = CGAL::Polygon_mesh_processing::does_self_intersect(m,
    CGAL::Polygon_mesh_processing::parameters::vertex_index_map(get(CGAL::vertex_point, m)));
  
  assert(intersecting_1 == intersecting_2);

  std::cerr << "does_self_intersect test took " << timer.time() << " sec." << std::endl;
  std::cerr
    << (intersecting_2 ? "There is a self-intersection." : "There are no self-intersections.")
    << std::endl;

  return 0;
}
