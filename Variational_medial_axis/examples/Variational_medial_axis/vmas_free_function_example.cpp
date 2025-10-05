#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/extract_variational_medial_skeleton.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;
using Medial_Skeleton = CGAL::Medial_Skeleton<Mesh>;

int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  Medial_Skeleton skeleton = extract_variational_medial_skeleton(
      mesh, CGAL::parameters::number_of_iterations(1000)// number of max iterations
      .number_of_spheres(200) // target number of spheres
      .number_of_samples(20000) // number of surface samples
      .lambda(0.2) // lambda parameter for the optimization
      .random_seed(10) // random seed for point sampling
      .concurrency_tag(CGAL::Parallel_tag{}) // use parallel execution
      .acceleration_structure(CGAL::BVH_tag{}) // use BVH for acceleration
      .verbose(true)); // enable verbose output
  // Write skeleton to PLY file
  std::string output_filename = "skeleton.ply";
  CGAL::IO::write_PLY(skeleton, output_filename, CGAL::parameters::stream_precision(9));

}
