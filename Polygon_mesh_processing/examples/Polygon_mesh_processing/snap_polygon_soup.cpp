#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Real_timer.h>

#include <boost/container/small_vector.hpp>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  const std::string filename = argc == 1 ? CGAL::data_file_path("meshes/elephant.off")
                                         : std::string(argv[1]);

  const int grid_size = argc <= 2 ? 23
                                  : std::stoi(std::string(argv[2]));

  const std::string out_file = "rounded_soup.off";

  std::vector<typename Kernel::Point_3> points;
  std::vector<boost::container::small_vector<std::size_t, 3>> triangles;

  std::cout << "Snap rounding on " << filename << "\n";
  if (!CGAL::IO::read_polygon_soup(filename, points, triangles))
  {
    std::cerr << "Cannot read " << filename << "\n";
    return 1;
  }

  PMP::repair_polygon_soup(points, triangles);
  PMP::triangulate_polygons(points, triangles);

  std::cout << "#points = " << points.size() << " and #triangles = " << triangles.size() << std::endl;
  std::cout << "Is 2-manifold: " << PMP::is_polygon_soup_a_polygon_mesh(triangles) << std::endl;

  std::vector<std::pair<std::size_t, std::size_t>> pairs_of_intersecting_triangles;
  PMP::triangle_soup_self_intersections(points, triangles, std::back_inserter(pairs_of_intersecting_triangles));
  std::cout << "Nb of pairs of intersecting triangles: " << pairs_of_intersecting_triangles.size() << std::endl;

  CGAL::Real_timer t;
  t.start();
  bool success=PMP::autorefine_triangle_soup(points, triangles, CGAL::parameters::apply_iterative_snap_rounding(true).erase_all_duplicates(false).concurrency_tag(CGAL::Parallel_if_available_tag()).snap_grid_size(grid_size).number_of_iterations(15));
  t.stop();

  std::cout << "\nOutput:" << std::endl;
  std::cout << "#points = " << points.size() << " and #triangles = " << triangles.size() << " in " << t.time() << " sec." << std::endl;
  if(success)
    std::cout << "Does self-intersect: " << PMP::does_triangle_soup_self_intersect<CGAL::Parallel_if_available_tag>(points, triangles) << std::endl;
  else
    std::cout << "ROUNDING FAILED" << std::endl;

  CGAL::IO::write_polygon_soup(out_file, points, triangles, CGAL::parameters::stream_precision(17));
  std::cout << "Is 2-manifold: " << PMP::orient_polygon_soup(points, triangles) << "\n\n"  << std::endl;
  return 0;
}
