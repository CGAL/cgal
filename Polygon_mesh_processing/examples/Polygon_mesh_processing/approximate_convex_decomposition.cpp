#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/approximate_convex_decomposition.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Point = K::Point_3;

using Convex_hull = std::pair<std::vector<Point>, std::vector<std::array<unsigned int, 3> > >;
using Mesh = CGAL::Surface_mesh<Point>;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/knot2.off");
  std::cout << filename << std::endl;

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::vector<Convex_hull> convex_hulls;
  PMP::approximate_convex_decomposition(mesh, std::back_inserter(convex_hulls),
    CGAL::parameters::maximum_depth(10)
    .volume_error(0.1)
    .maximum_number_of_convex_hulls(9)
    .split_at_concavity(true)
    .maximum_number_of_voxels(1000000)
    .concurrency_tag(CGAL::Parallel_if_available_tag()));

  for (std::size_t i = 0;i<convex_hulls.size();i++) {
    const Convex_hull& ch = convex_hulls[i];
    CGAL::IO::write_polygon_soup(std::to_string(i) + ".off", ch.first, ch.second);
  }

  return 0;
}
