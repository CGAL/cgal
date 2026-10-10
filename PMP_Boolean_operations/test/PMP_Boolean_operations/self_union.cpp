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
using Triangle = boost::container::small_vector<std::size_t, 3>;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

int main(int argc, char** argv)
{
  std::vector<Exact_point> points;
  std::vector<Triangle> triangles;

  const std::string filename = argc == 1 ? CGAL::data_file_path("meshes/elephant.off")
                                         : std::string(argv[1]);

  if (!CGAL::IO::read_polygon_soup(filename, points, triangles))
  {
    std::cerr << "Error reading " << filename << "\n";
    return 1;
  }

  CGAL::Real_timer timer;
  timer.start();
  PMP::compute_self_union(points, triangles, CGAL::parameters::concurrency_tag(Concurrency_tag()));
  timer.stop();
  std::cout << "Done in " << timer.time() << "\n";

  CGAL::IO::write_OFF("self_union.off", points, triangles, params::stream_precision(17));

  return 0;
}
