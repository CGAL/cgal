#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/pig.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  bool intersecting = PMP::does_self_intersect(mesh,
      PMP::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)));

  std::cout
    << (intersecting ? "There are self-intersections." : "There is no self-intersection.")
    << std::endl;

  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  PMP::self_intersections(mesh,
    std::back_inserter(intersected_tris),
    PMP::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)));

  std::cout << intersected_tris.size() << " pairs of triangles intersect." << std::endl;
  
  return 0;
}
