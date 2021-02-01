#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_2.h>
#include <CGAL/Projection_on_sphere_traits_3.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt  K;

// The kernel below cannot represent perfectly all points of the sphere
// and thus a separation mecanism is needed to ensure that no points are hidden

// typedef CGAL::Exact_predicates_inexact_constructions_kernel          K;

typedef K::FT                                                        FT;
typedef K::Point_3                                                   Point;

typedef CGAL::Projection_on_sphere_traits_3<K>                       Traits;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Traits>             DToS2;

int main(int argc, char** argv)
{
  std::cout.precision(17);

  const char* filename = (argc > 1) ? argv[1] : "data/poste_france.data";

  int n;
  std::vector<Point> points;
  double x, y, z;

  std::ifstream in(filename);
  in >> n;
  while(in >> x >> y >> z)
    points.emplace_back(x, y, z);

  // Add an extra point that would be too close to 'p' with a basic kernel such as CGAL::EPICK,
  const Point& p = points.back();
  const FT tiny = std::numeric_limits<double>::min();
  Point p_bis(p.x() + tiny, p.y() - tiny, p.z() + tiny);

  points.push_back(p_bis);

  std::cout << points.size() << " points in input" << std::endl;

  Traits traits(Point(0, 0, 0), 100);
  DToS2 dtos(points.begin(), points.end(), traits);

  std::cout << dtos.number_of_vertices() << " nv" << std::endl;
  std::cout << dtos.number_of_faces() << " nf" << std::endl;
  std::cout << dtos.number_of_ghost_faces() << " gf" << std::endl;

  std::ofstream out("result.off");
  out.precision(17);
  CGAL::write_OFF(out, dtos);
}
