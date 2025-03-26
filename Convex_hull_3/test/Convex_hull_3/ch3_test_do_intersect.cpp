#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Convex_hull_3/predicates.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3                                               Point_3;
typedef CGAL::Surface_mesh<Point_3>                              Mesh;

void test_cube()
{
  std::vector<Point_3> cube;
  for(int x=0; x<2; ++x)
    for(int y=0; y<2; ++y)
      for(int z=0; z<2; ++z)
          cube.push_back(Point_3(x,y,z));

  std::vector<Point_3> inside(1, Point_3(0.25,0.25,0.20));
  std::vector<Point_3> outside(1, Point_3(-0.25,0.25,0.25));

  bool res1=CGAL::Convex_hull_3::do_intersect(cube, inside);
  std::cout << "Do intersect inside: " << res1 << std::endl;
  bool res2=CGAL::Convex_hull_3::do_intersect(cube, outside);
  std::cout << "Do intersect outside: " << res2 << std::endl;
}

void test_half_sphere()
{
  std::vector<Point_3> half_sphere;
  // constexpr K::FT eps(std::pow(2,-40));
  constexpr double pi=3.14159265358979323846;
  for(double phi=25./16.; phi>0; phi-=1./4.)
    for(double theta=0; theta<2*pi; theta+=0.25)
          half_sphere.push_back(Point_3(std::sin(phi) * std::cos(theta),
                                        std::sin(phi) * std::sin(theta),
                                        std::cos(phi)));

  std::vector<Point_3> outside(1, Point_3(0.4,0.4,std::nextafter(std::cos(25./16),0)));
  std::vector<Point_3> inside(1, Point_3(0.4,0.4,std::nextafter(std::cos(25./16),1)));

  //print the number of extreme vertices
  std::cout << "There are " << half_sphere.size() << " extreme vertices in the input mesh.\n";
  // std::cout << bbox_3(half_sphere.begin(), half_sphere.end()) << std::endl;

  bool res1=CGAL::Convex_hull_3::do_intersect(half_sphere, inside);
  std::cout << "Do intersect inside: " << res1 << std::endl;
  bool res2=CGAL::Convex_hull_3::do_intersect(half_sphere, outside);
  std::cout << "Do intersect outside: " << res2 << std::endl;
}

int main()
{
  std::cout << std::setprecision(17);
  test_cube();
  test_half_sphere();

  return 0;
}
