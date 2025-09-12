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

  CGAL::Variational_medial_axis<Mesh, CGAL::Parallel_if_available_tag> vmas(mesh);
  // Compute medial skeleton with custom parameters
  vmas.compute_variational_medial_axis_sampling(
      CGAL::parameters::number_of_iterations(1000) // number of max iterations
      .number_of_spheres(200)// target number of spheres
      .lambda(0.2) // lambda parameter for the optimization
      .verbose(true)); // enable verbose output

  // add additional 3 sphere, this function will update the medial skeleton automatically
  vmas.add_spheres(3);
  // update the medial skeleton with 20 iterations
  vmas.update(20);

  std::cout << "Exporting skeleton..." << std::endl;
  auto skeleton = vmas.export_skeleton();

  // Write skeleton to PLY file
  std::string output_filename = "skeleton.ply";
  CGAL::IO::write_PLY(skeleton, output_filename, CGAL::parameters::stream_precision(9));

}

