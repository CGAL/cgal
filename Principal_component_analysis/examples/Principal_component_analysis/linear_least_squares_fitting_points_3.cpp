// Example program for the linear_least_square_fitting function on set of points in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line_3;
typedef K::Point_3           Point_3;

int main()
{
  std::list<Point_3> points;
  points.push_back(Point_3(1.0,0.0,0.0));
  points.push_back(Point_3(2.0,0.0,0.0));
  points.push_back(Point_3(3.0,0.0,0.0));

  Line_3 line;
  linear_least_squares_fitting_3(points.begin(),points.end(),line,CGAL::PCA_dimension_0_tag());

  return 0;
}
