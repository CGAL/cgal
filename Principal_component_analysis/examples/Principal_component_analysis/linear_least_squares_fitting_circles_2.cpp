// Example program for linear least squares fitting of 2D circles and disks
#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_2            Line;
typedef K::Point_2           Point;
typedef K::Circle_2          Circle;

int main()
{
	// generate a set of 2D circles
  std::list<Circle> circles;
  circles.push_back(Circle(Point(1.0,2.0),16.0));
  circles.push_back(Circle(Point(3.0,4.0),25.0));

	// fit line to circles
  Line line;
  linear_least_squares_fitting_2(circles.begin(),circles.end(),line,CGAL::Dimension_tag<1>());

	// fit line to disks
  linear_least_squares_fitting_2(circles.begin(),circles.end(),line,CGAL::Dimension_tag<2>());

  return 0;
}
