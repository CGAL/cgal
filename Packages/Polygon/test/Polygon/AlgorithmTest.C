#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <CGAL/std/fstream>
#include <CGAL/std/vector>
#include <CGAL/std/cstdlib>

//------------------------------------------------------------------------------//
//      test the polygon algorithms for a specific choice of the
//      representation class R and the point class Point
//------------------------------------------------------------------------------//
template <class R, class Point>
void test_polygon(const R&, const Point&, const char* FileName)
{
  typedef typename CGAL_STD::vector<Point>::iterator iterator;

  // read a point and a polygon from the file 'FileName'
  CGAL_STD::ifstream from(FileName);
  if (!from) {
    CGAL_STD::cerr << "could not open file " << FileName << "!" << endl;
    exit(1);
  }
  CGAL_set_ascii_mode(from);

  Point point;
  CGAL_STD::vector<Point> polygon;

  from >> point;
  copy(istream_iterator<Point, ptrdiff_t>(from),
       istream_iterator<Point, ptrdiff_t>(),
       back_inserter(polygon)
  );

    // show the point and the polygon
    CGAL_STD::cout << "point: " << point << endl;
    CGAL_STD::cout << "polygon:" << endl;
    copy(polygon.begin(), polygon.end(),
	ostream_iterator<Point>(CGAL_STD::cout, "\n"));
    CGAL_STD::cout << endl;

    iterator left =
	CGAL_left_vertex_2(polygon.begin(), polygon.end());
    iterator right =
	CGAL_right_vertex_2(polygon.begin(), polygon.end());
    iterator top =
	CGAL_top_vertex_2(polygon.begin(), polygon.end());
    iterator bottom =
	CGAL_bottom_vertex_2(polygon.begin(), polygon.end());
    bool simple =
	CGAL_is_simple_2(polygon.begin(), polygon.end());
    bool convex =
	CGAL_is_convex_2(polygon.begin(), polygon.end());
    CGAL_Bbox_2 bbox =
	CGAL_bbox_2(polygon.begin(), polygon.end());
    typename R::FT area = 0;
    CGAL_area_2(polygon.begin(), polygon.end(), area);
    CGAL_Bounded_side bside =
	 CGAL_bounded_side_2(polygon.begin(), polygon.end(), point);
    CGAL_Oriented_side oside =
	CGAL_oriented_side_2(polygon.begin(), polygon.end(), point);
    CGAL_Orientation orientation =
	CGAL_orientation_2(polygon.begin(), polygon.end());

    CGAL_STD::cout << "left   = " << *left << endl;
    CGAL_STD::cout << "right  = " << *right << endl;
    CGAL_STD::cout << "top    = " << *top << endl;
    CGAL_STD::cout << "bottom = " << *bottom << endl;
    CGAL_STD::cout << "the polygon is ";
    if (!simple) CGAL_STD::cout << "not ";
    CGAL_STD::cout << "simple" << endl;
    CGAL_STD::cout << "the polygon is ";
    if (!convex) CGAL_STD::cout << "not ";
    CGAL_STD::cout << "convex" << endl;
    CGAL_STD::cout << "the bounding box is " << bbox << endl;
    CGAL_STD::cout << "the area is " << area << endl;

    switch (bside) {
    case CGAL_ON_BOUNDED_SIDE:
	CGAL_STD::cout << "the point is on bounded side" << endl;
	break;
    case CGAL_ON_BOUNDARY:
	CGAL_STD::cout << "the point is on the boundary" << endl;
	break;
    case CGAL_ON_UNBOUNDED_SIDE:
	CGAL_STD::cout << "the point is on the unbounded side" << endl;
	break;
    }

    switch (oside) {
    case CGAL_ON_NEGATIVE_SIDE:
	CGAL_STD::cout << "the point is on the negative side" << endl;
	break;
    case CGAL_ON_ORIENTED_BOUNDARY:
	CGAL_STD::cout << "the point is on the oriented boundary" << endl;
	break;
    case CGAL_ON_POSITIVE_SIDE:
	CGAL_STD::cout << "the point is on the positive side" << endl;
	break;
    }

    switch(orientation) {
    case CGAL_CLOCKWISE:
	CGAL_STD::cout << "the orientation is clockwise" << endl;
	break;
    case CGAL_COUNTERCLOCKWISE:
	CGAL_STD::cout << "the orientation is counter clockwise" << endl;
	break;
    case CGAL_COLLINEAR:
	CGAL_STD::cout << "the orientation is collinear" << endl;
	break;
    }
}

//-----------------------------------------------------------------------//
//                          main
//-----------------------------------------------------------------------//

int main()
{
  CGAL_set_pretty_mode(CGAL_STD::cout);

  CGAL_STD::cout << endl;
  CGAL_STD::cout << "--------------------------------------------------------" << endl;
  CGAL_STD::cout << "-   testing polygon algorithms with cartesian/double   -" << endl;
  CGAL_STD::cout << "--------------------------------------------------------" << endl;
  CGAL_STD::cout << endl;
  typedef CGAL_Cartesian<double> R1;
//  CGAL_Polygon_traits_2<R1> traits1;

  typedef CGAL_Point_2<R1> Point1;
  test_polygon(R1(), Point1(), "data/polygon_cartesian.dat");

  CGAL_STD::cout << endl;
  CGAL_STD::cout << "--------------------------------------------------------" << endl;
  CGAL_STD::cout << "-   testing polygon algorithms with homogeneous/double -" << endl;
  CGAL_STD::cout << "--------------------------------------------------------" << endl;
  CGAL_STD::cout << endl;
  typedef CGAL_Homogeneous<double> R2;
  typedef CGAL_Point_2<R2> Point2;
  test_polygon(R2(), Point2(), "data/polygon_homogeneous.dat");

  return 0;
}

