// Example program for the linear_least_square_fitting function on a set of rectangles in 2D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_2            Line;
typedef K::Point_2           Point;
typedef K::Iso_rectangle_2   Iso_rectangle;

int main()
{
  std::list<Iso_rectangle> rectangles;
  rectangles.push_back(Iso_rectangle(Point(0.0,0.0),Point(4.0,8.0)));
  rectangles.push_back(Iso_rectangle(Point(4.0,8.0),Point(0.0,16.0)));

  Line line;

  // fit whole rectangles
  linear_least_squares_fitting_2(rectangles.begin(),rectangles.end(),line,CGAL::Dimension_tag<2>());
  
  // fit rectangle edges
  linear_least_squares_fitting_2(rectangles.begin(),rectangles.end(),line,CGAL::Dimension_tag<1>());
  
  // fit rectangle vertices
  linear_least_squares_fitting_2(rectangles.begin(),rectangles.end(),line,CGAL::Dimension_tag<0>());

  return 0;
}
