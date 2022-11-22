#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_ratio_stop_predicate.h>

#include <chrono>
#include <fstream>
#include <iostream>

typedef CGAL::Simple_cartesian<double>               Kernel;
typedef Kernel::Point_3                              Point_3;
typedef CGAL::Surface_mesh<Point_3>                  Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;

int main(int argc, char** argv)
{
  Surface_mesh surface_mesh;
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube-meshed.off");
  std::ifstream is(filename);
  if(!is || !(is >> surface_mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::is_triangle_mesh(surface_mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // In this example, the simplification stops when the number of undirected edges
  // drops below 10% of the initial count
  double stop_ratio = (argc > 2) ? std::stod(argv[2]) : 0.1;
  SMS::Edge_count_ratio_stop_predicate<Surface_mesh> stop(stop_ratio);

  int r = SMS::edge_collapse(surface_mesh, stop);

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

  std::cout << "\nFinished!\n" << r << " edges removed.\n" << surface_mesh.number_of_edges() << " final edges.\n";
  std::cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count() << "ms" << std::endl;

  CGAL::IO::write_polygon_mesh((argc > 3) ? argv[3] : "out.off", surface_mesh, CGAL::parameters::stream_precision(17));

  return EXIT_SUCCESS;
}
