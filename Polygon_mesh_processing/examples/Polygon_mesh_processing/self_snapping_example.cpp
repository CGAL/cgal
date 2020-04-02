#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel           Kernel;

typedef Kernel::FT                                                    FT;
typedef Kernel::Point_3                                               Point_3;
typedef Kernel::Vector_3                                              Vector_3;
typedef CGAL::Surface_mesh<Point_3>                                   Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor          vertex_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  std::cout.precision(17);

  if(argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " input_mesh tolerance" << std::endl;
    return EXIT_FAILURE;
  }

  Surface_mesh sm;
  std::ifstream in(argv[1]);
  if(!in || in >> sm)
  {
    std::cerr << "Problem loading the input data" << std::endl;
    return EXIT_FAILURE;
  }

  Surface_mesh::Property_map<vertex_descriptor, double> tolerance_map;
  tolerance_map = sm.add_property_map<vertex_descriptor, double>("v:t").first;
  for(vertex_descriptor v : vertices(sm))
    put(tolerance_map, v, std::atof(argv[3]));

  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  // Snap
  std::size_t nb_snapped = PMP::experimental::snap_borders<CGAL::Parallel_tag>(sm, tolerance_map);
  std::cout << "#snapped: " << nb_snapped << std::endl;

  std::chrono::steady_clock::time_point snap_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (snap): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(snap_time - start_time).count()
            << "ms" << std::endl;

  // Stitch
  std::cout << "Stitch, #ne: " << edges(sm).size() << std::endl;
  PMP::stitch_borders(sm);

  std::chrono::steady_clock::time_point stitch_time = std::chrono::steady_clock::now();
  std::cout << "Time elapsed (stitch): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(stitch_time - snap_time).count()
            << "ms" << std::endl;

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
  std::cout << "Total time elapsed: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
            << "ms" << std::endl;

  std::cout << "#border: " << PMP::number_of_borders(sm) << std::endl;
  std::cout << "Done!" << std::endl;

  std::ofstream("snapped.off") << std::setprecision(17) << sm;

  return EXIT_SUCCESS;
}
