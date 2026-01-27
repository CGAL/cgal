#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/tangential_relaxation.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;


namespace PMP = CGAL::Polygon_mesh_processing;


int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/pig.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  unsigned int nb_iter = (argc > 2) ? std::stoi(std::string(argv[2])) : 10;

  std::cout << "Relax...";

  PMP::tangential_relaxation(mesh, CGAL::parameters::number_of_iterations(nb_iter));

  std::cout << "done." << std::endl;

  return 0;
}
