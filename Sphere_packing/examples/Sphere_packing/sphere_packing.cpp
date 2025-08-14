#include <CGAL/Sphere_packing.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, const char **argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/bear.off");

  Mesh mesh;
  if (!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }
  return 0;
}
