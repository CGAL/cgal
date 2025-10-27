#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>

#include <CGAL/Surface_mesh.h>

#include <boost/container/small_vector.hpp>
#include <CGAL/Polygon_mesh_processing/orientation.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
typedef typename Kernel::Point_3 Point_3;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv)
{
  if(argc<4){
    std::cout << "Invalid argument" << std::endl;
    return 1;
  }
  const std::string filename = std::string(argv[1]);
  const int grid_size = std::stoi(std::string(argv[2]));
  const bool erase_duplicate = std::stoi(argv[3])==1;

  std::vector<Point_3> points;
  std::vector<boost::container::small_vector<std::size_t, 3>> triangles;

  CGAL::Bbox_3 bb = CGAL::bbox_3(points.begin(), points.end());
  double diag_length=std::sqrt((bb.xmax()-bb.xmin())*(bb.xmax()-bb.xmin()) + (bb.ymax()-bb.ymin())*(bb.ymax()-bb.ymin()) + (bb.zmax()-bb.zmin())*(bb.zmax()-bb.zmin()));
  if (!CGAL::IO::read_polygon_soup(filename, points, triangles))
  {
    std::cerr << "Cannot read " << filename << "\n";
    return 1;
  }

  std::vector<Point_3> input_points(points.begin(), points.end());

  PMP::autorefine_triangle_soup(points, triangles, CGAL::parameters::apply_iterative_snap_rounding(true).erase_all_duplicates(erase_duplicate).concurrency_tag(CGAL::Parallel_if_available_tag()).snap_grid_size(grid_size).number_of_iterations(15));


  std::cout << "{" <<
    "\"Nb_output_points\": \"" << points.size() << "\",\n" <<
    "\"Nb_output_triangles\": \"" << triangles.size() << "\",\n" <<
    "\"Is_2_manifold\": \"" << (PMP::orient_polygon_soup(points, triangles)?"True":"False") << "\",\n";
  CGAL::Surface_mesh<Point_3> sm;
  PMP::polygon_soup_to_polygon_mesh(points, triangles, sm);

  std::cout << std::setprecision(17) <<
   "\"Hausdorff_distance_output_to_input_(divide_by_bbox_diag)\": \"" << PMP::max_distance_to_triangle_mesh<CGAL::Parallel_if_available_tag>(input_points, sm) / diag_length << "\",\n" <<
   "\"Closed_output\": \"" << (CGAL::is_closed(sm)?"True":"False") << "\",\n" <<
   "\"Ouput_bound_a_volume\": \"" << (PMP::does_bound_a_volume(sm)?"True":"False") << "\"\n}"
  << std::endl;

  return 0;
}