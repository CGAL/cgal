#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/variational_medial_axis_sampling.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;

int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  CGAL::variational_medial_axis_sampling(mesh);

}

