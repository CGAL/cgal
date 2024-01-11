#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point;

typedef CGAL::Surface_mesh<Point>                           Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  Mesh mesh;

  const std::string filename = argc == 1 ? CGAL::data_file_path("meshes/elephant.off")
                                         : std::string(argv[1]);
  CGAL::IO::read_polygon_mesh(filename, mesh);

  PMP::autorefine(mesh);

  CGAL::IO::write_polygon_mesh("autorefined.off", mesh, CGAL::parameters::stream_precision(17));

  return 0;
}
