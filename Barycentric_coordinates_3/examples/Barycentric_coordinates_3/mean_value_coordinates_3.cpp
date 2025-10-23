#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/generators.h>
#include <CGAL/Barycentric_coordinates_3/Mean_value_coordinates_3.h>

using Kernel =  CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = Kernel::FT;
using Point_3 =  Kernel::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;
using CP3 = CGAL::Barycentric_coordinates::Computation_policy_3;

int main() {

  Surface_mesh concave;

  const Point_3 p0(0, 3, 0);
  const Point_3 p1(1, 1, 0);
  const Point_3 p2(3, 0, 0);
  const Point_3 p3(0, 0, 0);
  const Point_3 p4(0, 0, 3);
  const Point_3 p5(0, 3, 3);
  const Point_3 p6(1, 1, 3);
  const Point_3 p7(3, 0, 3);

  CGAL::make_hexahedron(p0, p1, p2, p3, p4, p5, p6, p7, concave,
                        CGAL::parameters::do_not_triangulate_faces(false));

  std::vector<FT> coords;
  std::vector<Point_3> queries{
    Point_3(FT(1)/FT(2), FT(1)/FT(2), FT(1)), Point_3(FT(1)/FT(3), FT(1)/FT(3), FT(2)), // Only points in the kernel
    Point_3(FT(4)/FT(3), FT(1)/FT(3), FT(1)), Point_3(FT(4)/FT(3), FT(1)/FT(3), FT(2)),
    Point_3(FT(1)/FT(3), FT(4)/FT(3), FT(1)), Point_3(FT(1)/FT(3), FT(4)/FT(3), FT(2))};

  std::cout << std::endl << "Mean value coordinates : " << std::endl << std::endl;

  for (const auto& query : queries){

    coords.clear();
    CGAL::Barycentric_coordinates::mean_value_coordinates_3(
      concave, query, std::back_inserter(coords), CP3::FAST);

    // Output mean value coordinates.
    for (std::size_t i = 0; i < coords.size() -1; ++i) {
      std::cout << coords[i] << ", ";
    }
    std::cout << coords[coords.size() -1] << std::endl;
  }
  std::cout << std::endl;


  return EXIT_SUCCESS;
}
