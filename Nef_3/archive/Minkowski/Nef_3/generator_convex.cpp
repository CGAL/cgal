#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/point_generators_3.h>

typedef CGAL::Gmpz NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Random_points_in_cube_3<Point_3> Point_source;

int main(int argc, char* argv[]) {

  CGAL_assertion(argc>1 && argc < 5);

  int n = std::atoi(argv[1]);
  int s = (argc > 2 ? std::atoi(argv[2]) : 10000);
  unsigned int seed = (argc > 3 ? std::atoi(argv[3]) : 12345);

  CGAL::Random R(seed);
  Point_source PS(s, R);
  std::list<Point_3> points;

  for(int i=0; i<n; ++i)
    points.push_back(*PS++);

  Polyhedron P;
  convex_hull_3( points.begin(), points.end(), P);

  std::cout << P;
}
