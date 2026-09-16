#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/subdivision_method_3.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>          Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>     Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

int main(int argc, char* argv[]) {

  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/quad.off");
  const unsigned int iter = (argc > 2) ? std::stoi(argv[2]) : 3;

  Surface_mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  // Loop_subdivision can also be used for pure triangle meshes
  CGAL::Subdivision_method_3::CatmullClark_subdivision(mesh, params::number_of_iterations(iter)
                                                                    .do_not_modify_geometry(true));

  std::ofstream out("out.off");
  out << mesh;

  return 0;
}
