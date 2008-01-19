// Example program for linear least squares fitting of a 2D point set
#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_2            Line;
typedef K::Point_2           Point;

int main()
{
  std::list<Point> points;
  points.push_back(Point(1.0,2.0));
  points.push_back(Point(3.0,4.0));
  points.push_back(Point(5.0,6.0));

	// fit a line 
  Line line;
  linear_least_squares_fitting_2(points.begin(),points.end(),line,CGAL::PCA_dimension_0_tag());

  return 0;
}
