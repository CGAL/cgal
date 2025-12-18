#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <boost/container/small_vector.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef typename Kernel::Point_3 Point_3;
namespace PMP = CGAL::Polygon_mesh_processing;

enum EXIT_CODES { VALID_OUTPUT=0,
                  INVALID_INPUT=1,
                  ROUNDING_FAILED=2,
                  SELF_INTERSECTING_OUTPUT=3,
                  SIGSEGV=10,
                  SIGSABRT=11,
                  SIGFPE=12,
                  TIMEOUT=13
                  };

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

  if (!CGAL::IO::read_polygon_soup(filename, points, triangles) || points.size()==0 || triangles.size()==0)
  {
    return INVALID_INPUT;
  }

  PMP::repair_polygon_soup(points, triangles);
  PMP::triangulate_polygons(points, triangles);

  bool success=PMP::autorefine_triangle_soup(points, triangles, CGAL::parameters::apply_iterative_snap_rounding(true).erase_all_duplicates(erase_duplicate).concurrency_tag(CGAL::Parallel_if_available_tag()).snap_grid_size(grid_size).number_of_iterations(15));

  if(!success)
    return ROUNDING_FAILED;
  if( PMP::does_triangle_soup_self_intersect<CGAL::Parallel_if_available_tag>(points, triangles) )
    return SELF_INTERSECTING_OUTPUT;

  return VALID_OUTPUT;
}
