#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor    vertex_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  if(argc != 4)
  {
    std::cerr << "Usage: " << argv[0] << " movable_mesh fixed_mesh tolerance" << std::endl;
    return EXIT_FAILURE;
  }

  Surface_mesh movable_mesh, fixed_mesh;

  std::ifstream in_m(argv[1]);
  if(!in_m || in_m >> movable_mesh)
  {
    std::cerr << "Problem loading the input data" << std::endl;
    return EXIT_FAILURE;
  }

  std::ifstream in_f(argv[2]);
  if(!in_f || in_f >> fixed_mesh)
  {
    std::cerr << "Problem loading the input data" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Movable mesh: " << num_vertices(movable_mesh) << std::endl;
  std::cout << "Fixed mesh: " << num_vertices(fixed_mesh) << std::endl;

  const double tolerance = std::atof(argv[3]);

  Surface_mesh::Property_map<vertex_descriptor, double> movable_tolerance_map;
  movable_tolerance_map = movable_mesh.add_property_map<vertex_descriptor, double>("v:t").first;
  for(vertex_descriptor v : vertices(movable_mesh))
    put(movable_tolerance_map, v, tolerance);

  Surface_mesh::Property_map<vertex_descriptor, double> fixed_tolerance_map;
  fixed_tolerance_map = fixed_mesh.add_property_map<vertex_descriptor, double>("v:t").first;
  for(vertex_descriptor v : vertices(fixed_mesh))
    put(fixed_tolerance_map, v, tolerance);

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // Choice of named parameters indicates that:
  // - We want to simplify the boundary of the first mesh
  // - The geometry of the second mesh cannot change
  std::size_t nb_snapped = PMP::experimental::snap_borders(movable_mesh, movable_tolerance_map,
                                                           fixed_mesh, fixed_tolerance_map,
                                                           CGAL::parameters::do_simplify_border(true),
                                                           CGAL::parameters::do_lock_mesh(true));

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

  std::cout << "Time elapsed: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
            << "ms" << std::endl;

  std::cout << "#Snapped: " << nb_snapped << std::endl;

  std::ofstream("snapped_movable.off") << std::setprecision(17) << movable_mesh;
  std::ofstream("snapped_fixed.off") << std::setprecision(17) << fixed_mesh;

  return EXIT_SUCCESS;
}
