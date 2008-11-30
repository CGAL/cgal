// Example program for linear least squares fitting of 3D segments
#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line;
typedef K::Plane_3           Plane;
typedef K::Point_3           Point;
typedef K::Segment_3         Segment;

int main()
{
	Point a(1.0,2.0,3.0);
	Point b(4.0,5.0,6.0);
	Point c(7.0,8.0,9.0);
  std::list<Segment> segments;
  segments.push_back(Segment(a,b));
  segments.push_back(Segment(a,c));

  Line line;
  Plane plane;

	// fit line and plane to whole segments
  linear_least_squares_fitting_3(segments.begin(),segments.end(),line, CGAL::Dimension_tag<1>());
  linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,CGAL::Dimension_tag<1>());

	// fit line and plane to segment end points
  linear_least_squares_fitting_3(segments.begin(),segments.end(),line, CGAL::Dimension_tag<0>());
  linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,CGAL::Dimension_tag<0>());

  return 0;
}
