// Example program for linear least squares fitting of 2D segments
#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_2            Line;
typedef K::Point_2           Point;
typedef K::Segment_2         Segment;

int main()
{
	Point a(1.0,2.0);
	Point b(3.0,4.0);
	Point c(5.0,6.0);
  std::list<Segment> segments;
  segments.push_back(Segment(a,b));
  segments.push_back(Segment(a,c));

  Line line;
  
	// fit a line
  linear_least_squares_fitting_2(segments.begin(),segments.end(),line,CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_2(segments.begin(),segments.end(),line,CGAL::Dimension_tag<0>());

  return 0;
}
