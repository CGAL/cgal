// Example program for the linear_least_square_fitting function
// on a set of 3D triangles
#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <vector>

typedef double                      FT;
typedef CGAL::Simple_cartesian<FT>  K;
typedef K::Line_3                   Line;
typedef K::Plane_3                  Plane;
typedef K::Point_3                  Point;
typedef K::Triangle_3               Triangle;

int main(void)
{
  std::vector<Triangle> triangles;
  Point a(1.0,2.0,3.0);
  Point b(4.0,0.0,6.0);
  Point c(7.0,8.0,9.0);
  Point d(8.0,7.0,6.0);
  Point e(5.0,3.0,4.0);
  triangles.push_back(Triangle(a,b,c));
  triangles.push_back(Triangle(a,b,d));
  triangles.push_back(Triangle(d,e,c));

  Line line;
  Plane plane;

  // fit plane to whole triangles
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<2>());
  
  // fit line to triangle vertices
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::Dimension_tag<0>());
  
  return 0;
}
