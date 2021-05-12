// Example program for linear least squares fitting of a 2D point set
#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <vector>

typedef double                      FT;
typedef CGAL::Simple_cartesian<FT>  K;
typedef K::Line_2                   Line;
typedef K::Point_2                  Point;

int main()
{
  std::vector<Point> points;
  points.push_back(Point(1.0,2.0));
  points.push_back(Point(3.0,4.0));
  points.push_back(Point(5.0,6.0));

  // fit line
  Line line;
  linear_least_squares_fitting_2(points.begin(),points.end(),line,CGAL::Dimension_tag<0>());

  return 0;
}
