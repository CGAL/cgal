// Example program for the linear_least_square_fitting function on set of segments in 3D

#include <CGAL/Cartesian.h>
#include <CGAL/linear_least_squares_fitting_3.h>

#include <list>

typedef double               FT;
typedef CGAL::Cartesian<FT>  K;
typedef K::Line_3            Line_3;
typedef K::Plane_3           Plane_3;
typedef K::Point_3           Point_3;
typedef K::Segment_3         Segment_3;

int main()
{
  std::list<Segment_3> segments;
  segments.push_back(Segment_3(Point_3(1.0,1.0,1.0),Point_3(2.0,2.0,2.0)));
  segments.push_back(Segment_3(Point_3(3.0,3.0,3.0),Point_3(8.0,8.0,8.0)));

  Line_3 line;
  Point_3 c;
  linear_least_squares_fitting_3(segments.begin(),segments.end(),line,c,CGAL::PCA_dimension_1_tag());

  Plane_3 plane;
  linear_least_squares_fitting_3(segments.begin(),segments.end(),plane,c,CGAL::PCA_dimension_1_tag());

  return 0;
}
