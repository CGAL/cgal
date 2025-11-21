#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Barycentric_coordinates_3/Wachspress_coordinates_3.h>

using Kernel =  CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = Kernel::FT;
using Point_3 =  Kernel::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;
using CP3 = CGAL::Barycentric_coordinates::Computation_policy_3;
using WP = CGAL::Barycentric_coordinates::Wachspress_coordinates_3<Surface_mesh, Kernel>;

int main() {

  CGAL::Random_points_in_sphere_3<Point_3> gen(1.0);
  std::vector<Point_3> points;

  const std::size_t number_of_points = 250;
  std::copy_n(gen, number_of_points, std::back_inserter(points));

  Surface_mesh sm;
  CGAL::convex_hull_3(points.begin(), points.end(), sm);
  const std::size_t number_of_vertices = num_vertices(sm);

  WP wp(sm, CP3::FAST_WITH_EDGE_CASES);

  std::cout << "Computed Wachspress coordinates: " << std::endl << std::endl;
  for(std::size_t i = 0; i < number_of_points; i++){

      std::vector<FT> coordinates;
      coordinates.reserve(number_of_vertices);
      wp(points[i], std::back_inserter(coordinates));

      std::cout << "Point " << i + 1 << ": " << std::endl;
      for(std::size_t j = 0; j < number_of_vertices; j++)
        std::cout << "Coordinate " << j + 1 << " = " << coordinates[j] << "; " << std::endl;
      std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}
