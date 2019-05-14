#include <CGAL/internal/disable_deprecation_warnings_and_errors.h>

// Example program for the linear_least_square_fitting function on a set of circles in 2D

#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>

#include <list>

typedef double                      FT;
typedef CGAL::Simple_cartesian<FT>  K;
typedef K::Line_2                   Line;
typedef K::Point_2                  Point;
typedef K::Circle_2                 Circle;

int main()
{
  std::cout << "Test 2D linear least squares fitting of circles"  << std::endl;

  std::list<Circle> circles;
  circles.push_back(Circle(Point(0.0, 0.0),9.0));
  circles.push_back(Circle(Point(0.0,10.0),49.0));
  circles.push_back(Circle(Point(10.0,0.0),49.0));

  Line line;
  Point centroid;

  linear_least_squares_fitting_2(circles.begin(),circles.end(),line,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_2(circles.begin(),circles.end(),line,CGAL::Dimension_tag<1>());

  linear_least_squares_fitting_2(circles.begin(),circles.end(),line,centroid,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_2(circles.begin(),circles.end(),line,centroid,CGAL::Dimension_tag<1>());

  return 0;
}
