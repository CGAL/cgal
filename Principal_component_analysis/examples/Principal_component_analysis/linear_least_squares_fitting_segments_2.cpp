// Example program for the linear_least_square_fitting function on set of segments in 2D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_2            Line_2;
typedef K::Point_2           Point_2;
typedef K::Segment_2         Segment_2;

int main()
{
  std::list<Segment_2> segments;
  segments.push_back(Segment_2(Point_2(0.0,1.0),Point_2(-1.0,0.0)));
  segments.push_back(Segment_2(Point_2(0.0,1.0),Point_2(1.0,0.0)));

  Line_2 line;
  Point_2 c;
  linear_least_squares_fitting_2(segments.begin(),segments.end(),line,c,CGAL::PCA_dimension_1_tag());

  return 0;
}
