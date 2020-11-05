#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>


#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                    Point;
typedef CGAL::Surface_mesh<Point>          Surface_mesh;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/P.off";
  const char* outfilename = (argc > 2) ? argv[2] : "P_tri.off";
  std::ifstream input(filename);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << "Not a valid off file." << std::endl;
    return 1;
  }

  CGAL::Polygon_mesh_processing::triangulate_faces(mesh);

  // Confirm that all faces are triangles.
  for(boost::graph_traits<Surface_mesh>::face_descriptor fit : faces(mesh))
    if (next(next(halfedge(fit, mesh), mesh), mesh)
        !=   prev(halfedge(fit, mesh), mesh))
      std::cerr << "Error: non-triangular face left in mesh." << std::endl;

  std::ofstream cube_off(outfilename);
  cube_off.precision(17);
  cube_off << mesh;

  return 0;
}
