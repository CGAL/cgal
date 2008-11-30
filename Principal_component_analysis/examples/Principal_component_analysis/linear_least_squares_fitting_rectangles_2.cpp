// Example program for the linear_least_square_fitting function on a set of rectangles in 2D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_2            Line_2;
typedef K::Point_2           Point_2;
typedef K::Iso_rectangle_2   Iso_rectangle_2;

int main()
{
  std::list<Iso_rectangle_2> Iso_rectangles;
  Iso_rectangles.push_back(Iso_rectangle_2(Point_2(0.0,0.0),Point_2(4.0,8.0)));
  Iso_rectangles.push_back(Iso_rectangle_2(Point_2(4.0,8.0),Point_2(0.0,16.0)));

  Line_2 line;
  Point_2 c;

  linear_least_squares_fitting_2(Iso_rectangles.begin(),Iso_rectangles.end(),line,c,CGAL::Dimension_tag<2>());
  linear_least_squares_fitting_2(Iso_rectangles.begin(),Iso_rectangles.end(),line,c,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_2(Iso_rectangles.begin(),Iso_rectangles.end(),line,c,CGAL::Dimension_tag<0>());

  return 0;
}
