#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_2.h>
#include <CGAL/Projection_on_sphere_traits_3.h>

#include <iostream>
#include <fstream>

template <typename Kernel>
void create_triangulation(const std::string& filename)
{
  typedef typename Kernel::FT                                          FT;

  typedef CGAL::Projection_on_sphere_traits_3<Kernel>                  Traits;
  typedef CGAL::Delaunay_triangulation_on_sphere_2<Traits>             DToS2;

  typedef typename Traits::Point_3                                     Point_3;

  std::cout << "\n-- Constructing triangulation with Kernel: " << typeid(Kernel).name() << " --" << std::endl;

  std::vector<Point_3> points;
  double x, y, z;

  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Invalid input file: " << filename << std::endl;
    return;
  }

  while(in >> x >> y >> z)
    points.emplace_back(x, y, z);

  // Add an extra point that would be too close to 'p' with a basic kernel such as CGAL::EPICK,
  const Point_3& p = points.back();
  const FT tiny = 100 * std::numeric_limits<double>::epsilon();
  points.emplace_back(p.x() + tiny, p.y() - tiny, p.z() + tiny);

  std::cout << "Adding point  " << points.back() << "\nvery close to " << p << std::endl;
  std::cout << "Squared distance between points " << CGAL::squared_distance(points.back(), p) << std::endl;
  std::cout << points.size() << " points in input" << std::endl;

  Traits traits(Point_3(0, 0, 0), 100); // centered on (0,0,0), with radius 100
  DToS2 dtos(points.begin(), points.end(), traits);

  std::cout << dtos.number_of_vertices() << " vertices" << std::endl;
  std::cout << dtos.number_of_faces() << " faces" << std::endl;
}

int main(int argc, char** argv)
{
  std::cout.precision(17);

  // This kernel CAN represent exactly all points of the sphere
  typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  EPECK_w_SQRT;

  // This kernel CANNOT represent exactly all points of the sphere
  // and thus a separation mechanism is needed to ensure that no points are hidden
  typedef CGAL::Exact_predicates_inexact_constructions_kernel          EPICK;

  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/poste_france.xyz");

  create_triangulation<EPICK>(filename);
  create_triangulation<EPECK_w_SQRT>(filename);

  return EXIT_SUCCESS;
}
