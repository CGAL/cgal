#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/variational_medial_axis_sampling.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;

void normalize_mesh(Mesh& mesh) {
  CGAL::Bbox_3 bbox;
  for(auto v : mesh.vertices())
    bbox = bbox + mesh.point(v).bbox();

  double cx = (bbox.xmin() + bbox.xmax()) / 2.0;
  double cy = (bbox.ymin() + bbox.ymax()) / 2.0;
  double cz = (bbox.zmin() + bbox.zmax()) / 2.0;
  double max_dim = std::max({bbox.xmax() - bbox.xmin(), bbox.ymax() - bbox.ymin(), bbox.zmax() - bbox.zmin()});

  for(auto v : mesh.vertices()) {
    K::Point_3 p = mesh.point(v);
    double nx = (p.x() - bbox.xmin()) / max_dim;
    double ny = (p.y() - bbox.ymin()) / max_dim;
    double nz = (p.z() - bbox.zmin()) / max_dim;
    mesh.point(v) = K::Point_3(nx, ny, nz);
  }
}

int main(int argc, char** argv)
{ 
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }
  normalize_mesh(mesh);
  
  CGAL::Variational_medial_axis<Mesh> vmas(mesh);
  vmas.init();
  // Compute medial axis with custom parameters
  //vmas.compute_variational_medial_axis_sampling(CGAL::parameters::number_of_iterations(2000).number_of_spheres(1000).lambda(1.0).concurrency_tag(CGAL::Parallel_tag{}));
  vmas.compute_variational_medial_axis_sampling(
      CGAL::parameters::number_of_iterations(1000).number_of_spheres(200).lambda(0.2).concurrency_tag(
      CGAL::Parallel_tag{}));
  // Export skeleton
  std::cout << "Exporting skeleton..." << std::endl;
  auto skeleton = vmas.export_skeleton();

  // Write skeleton to PLY file
  std::string output_filename = "skeleton.ply";
  skeleton.write_to_ply_file(output_filename);

}

