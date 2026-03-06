#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangle_soup_boolean_operations_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <boost/container/small_vector.hpp>

#include <fstream>

#include <CGAL/Real_timer.h>

#ifdef CGAL_LINKED_WITH_TBB
using Concurrency_tag = CGAL::Parallel_tag;
#else
using Concurrency_tag = CGAL::Sequential_tag;
#endif

typedef CGAL::Exact_predicates_exact_constructions_kernel    EPECK;
typedef EPECK::Point_3                                       Exact_point;
typedef CGAL::Surface_mesh<Exact_point>                      Exact_mesh;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;


void run(std::string f1, std::string f2)
{
  std::cout << "Running with " << f1 << " and " << f2 << "\n";
  using Triangle = boost::container::small_vector<std::size_t, 3>;

  std::vector<Exact_point> points_1, points_2;
  std::vector<Triangle> triangles_1, triangles_2;

  CGAL::IO::read_polygon_soup(f1, points_1, triangles_1);
  CGAL::IO::read_polygon_soup(f2, points_2, triangles_2);

{
  std::vector<Exact_point> points_res;
  std::vector<Triangle> triangles_res;
  CGAL::Real_timer timer;
  timer.start();
  PMP::compute_union<Concurrency_tag>(points_1, triangles_1, points_2, triangles_2, points_res, triangles_res);
  timer.stop();
  std::cout << "Done in " << timer.time() << "\n";

  CGAL::IO::write_OFF("union.off", points_res, triangles_res, params::stream_precision(17));
}
{
  std::vector<Exact_point> points_res;
  std::vector<Triangle> triangles_res;
  CGAL::Real_timer timer;
  timer.start();
  PMP::compute_intersection<Concurrency_tag>(points_1, triangles_1, points_2, triangles_2, points_res, triangles_res);
  timer.stop();
  std::cout << "Done in " << timer.time() << "\n";

  CGAL::IO::write_OFF("intersection.off", points_res, triangles_res, params::stream_precision(17));
}
{
  std::vector<Exact_point> points_res;
  std::vector<Triangle> triangles_res;
  CGAL::Real_timer timer;
  timer.start();
  PMP::compute_difference<Concurrency_tag>(points_1, triangles_1, points_2, triangles_2, points_res, triangles_res);
  timer.stop();
  std::cout << "Done in " << timer.time() << "\n";

  CGAL::IO::write_OFF("difference.off", points_res, triangles_res, params::stream_precision(17));
}
}

int main(int argc, char** argv)
{
  for (int i=0; i< (argc-1)/2; ++i)
    run(argv[2*i+1], argv[2*i+2]);
}
