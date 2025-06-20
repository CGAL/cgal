#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/variational_medial_axis_sampling.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;

int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/chair.off");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

 
  CGAL::Variational_medial_axis<Mesh, K> vmas(mesh);

  // Initialize with default parameters
  vmas.init();

  // Compute medial axis with custom parameters
  vmas.compute(CGAL::parameters::number_of_iterations(500).number_of_spheres(100).lambda(0.2));
 
  // Export skeleton
  std::cout << "Exporting skeleton..." << std::endl;
  auto skeleton = vmas.export_skeleton();

  // Write skeleton to PLY file
  std::string output_filename = "skeleton.ply";
  skeleton.write_to_ply_file(output_filename);

}

