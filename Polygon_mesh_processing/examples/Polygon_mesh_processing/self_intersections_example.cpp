#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/polygon_mesh_self_intersections.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

int main(int argc, char* argv[])
{
  char* filename = (argc > 1) ? argv[1] : "data/pig.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  bool intersecting =
    CGAL::Polygon_mesh_processing::is_self_intersecting<K>(mesh);

  std::cout
    << (intersecting ? "There are self-intersections." : "There is no self-intersection.")
    << std::endl;

  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  CGAL::Polygon_mesh_processing::self_intersections<K>
            (mesh, back_inserter(intersected_tris));

  std::cout << intersected_tris.size() << " pairs of triangles are intersecting." << std::endl;
  
  return 0;
}
