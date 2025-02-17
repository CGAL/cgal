#define  PMP_ROUNDING_VERTICES_IN_POLYGON_SOUP_VERBOSE
#define CGAL_PMP_AUTOREFINE_USE_DEFAULT_VERBOSE

#include <CGAL/Cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Real_timer.h>

#include <boost/container/small_vector.hpp>

#include <iostream>

#include <CGAL/Polygon_mesh_processing/internal/snap_polygon_soup.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Cartesian<double> Cartesian;
typedef Kernel::Point_3                                     Point;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const std::string filename = argc == 1 ? CGAL::data_file_path("meshes/elephant.off")
                                         : std::string(argv[1]);

  const int grid_size = argc == 2 ? 23
                                  : std::stoi(std::string(argv[2]));

  std::vector<Point> input_points;
  std::vector<boost::container::small_vector<std::size_t, 3>> input_triangles;
  if (!CGAL::IO::read_polygon_soup(filename, input_points, input_triangles))
  {
    std::cerr << "Cannot read " << filename << "\n";
    return 1;
  }
  PMP::repair_polygon_soup(input_points, input_triangles);
  PMP::triangulate_polygons(input_points, input_triangles);

  CGAL::Real_timer t;
  t.start();
  PMP::autorefine_triangle_soup(input_points, input_triangles, CGAL::parameters::do_snap(true).erase_all_duplicates(true).concurrency_tag(CGAL::Parallel_if_available_tag()).snap_grid_size(grid_size));
  t.stop();
  std::cout << "#points = " << input_points.size() << " and #triangles = " << input_triangles.size() << " in " << t.time() << " sec." << std::endl;
  
  std::vector<Cartesian::Point_3> output_points;
  std::cout << "Does self-intersect: " << PMP::does_triangle_soup_self_intersect(input_points, input_triangles) << std::endl;
  for(auto &p: input_points)
    output_points.emplace_back(CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));

  CGAL::IO::write_polygon_soup("rounded_soup.off", output_points, input_triangles, CGAL::parameters::stream_precision(17));

  return 0;
}
