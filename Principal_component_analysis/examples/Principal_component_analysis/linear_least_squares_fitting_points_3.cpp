// Example program for linear least squares fitting of 3D points
#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line;
typedef K::Plane_3           Plane;
typedef K::Point_3           Point;

int main()
{
  std::list<Point> points;
  points.push_back(Point(1.0,2.0,3.0));
  points.push_back(Point(4.0,5.0,6.0));
  points.push_back(Point(7.0,8.0,9.0));

	// fit line
  Line line;
  linear_least_squares_fitting_3(points.begin(),points.end(),line,CGAL::Dimension_tag<0>());

	// fit plane
  Plane plane;
  linear_least_squares_fitting_3(points.begin(),points.end(),plane,CGAL::Dimension_tag<0>());

  return 0;
}
