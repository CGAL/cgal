#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/approximate_convex_decomposition.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/convexity_check_3.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using EPICK = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;

namespace PMP = CGAL::Polygon_mesh_processing;

template<typename K>
void test(const std::string &filename, const std::vector<std::size_t> n) {
  using Point = typename K::Point_3;

  using Convex_hull = std::pair<std::vector<Point>, std::vector<std::array<unsigned int, 3> > >;
  using Mesh = CGAL::Surface_mesh<Point>;

  Mesh mesh;
  if (!PMP::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Invalid input." << std::endl;
    assert(false);
  }

  for (std::size_t i : n) {
    std::vector<Convex_hull> convex_hulls;
    PMP::approximate_convex_decomposition(mesh, std::back_inserter(convex_hulls),
      CGAL::parameters::maximum_depth(10)
      .volume_error(0.1)
      .maximum_number_of_convex_hulls(i)
      .split_at_concavity(true)
      .maximum_number_of_voxels(1000000)
      .concurrency_tag(CGAL::Parallel_if_available_tag()));
    std::cout << convex_hulls.size() << std::endl;
    assert(convex_hulls.size() == i);
    for (std::size_t j = 0; j < convex_hulls.size(); j++) {
      Mesh m;
      PMP::polygon_soup_to_polygon_mesh(convex_hulls[j].first, convex_hulls[j].second, m);

      assert(CGAL::is_strongly_convex_3(m, typename CGAL::Convex_hull_3::internal::Default_traits_for_Chull_3<Point>::type()));
    }
  }
}

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/knot2.off");
  std::cout << "Testing approximate convex decomposition with EPICK Kernel of " << filename << std::endl;
  test<EPICK>(filename, {3, 6, 9});
  std::cout << "Testing approximate convex decomposition with EPECK Kernel of " << filename << std::endl;
  test<EPECK>(filename, {5});

  std::cout << "done" << std::endl;

  return 0;
}
