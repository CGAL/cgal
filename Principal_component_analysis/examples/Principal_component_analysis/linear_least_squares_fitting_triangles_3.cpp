// Example program for the linear_least_square_fitting function
// on a set of 3D triangles
#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line;
typedef K::Plane_3           Plane;
typedef K::Point_3           Point;
typedef K::Triangle_3        Triangle;

int main(void)
{
  std::list<Triangle> triangles;
	Point a(1.0,2.0,3.0);
	Point b(4.0,5.0,6.0);
	Point c(7.0,8.0,9.0);
	Point d(8.0,7.0,6.0);
  triangles.push_back(Triangle(a,b,c));
  triangles.push_back(Triangle(a,b,d));

  Line line;
  Plane plane;

	// fit line and plane to whole triangles
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<2>());
  
	// fit line and plane to triangle edges
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<1>());
  
	// fit line and plane to triangle vertices
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),line, CGAL::Dimension_tag<0>());
  linear_least_squares_fitting_3(triangles.begin(),triangles.end(),plane,CGAL::Dimension_tag<0>());
  
  return 0;
}
