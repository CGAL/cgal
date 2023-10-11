#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/subdivision_method_3.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>          Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>     Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[]) {

  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/quad.off");

  Surface_mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  CGAL::Linear_mask_3<Surface_mesh> mask(&mesh);
  CGAL::Subdivision_method_3::PQQ(mesh, mask, CGAL::parameters::number_of_iterations(1));

  std::ofstream out("out.off");
  out << mesh;

  return 0;
}
