#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <fstream>
#include <vector>
#include <cstdlib>

using std::cout;
using std::endl;

//------------------------------------------------------------------------------//
//      test the polygon algorithms for a specific choice of the
//      representation class R and the point class Point
//------------------------------------------------------------------------------//
template <class R, class Point>
void test_polygon(const R&, const Point&, const char* FileName)
{
  typedef typename std::vector<Point>::iterator iterator;

  // read a point and a polygon from the file 'FileName'
  std::ifstream from(FileName);
  if (!from) {
    std::cerr << "could not open file " << FileName << "!" << endl;
    std::exit(1);
  }
  CGAL::set_ascii_mode(from);

  Point point;
  std::vector<Point> polygon;

  from >> point;
  std::copy(std::istream_iterator<Point>(from),
       std::istream_iterator<Point>(),
       std::back_inserter(polygon)
  );

    // show the point and the polygon
    cout << "point: " << point << endl;
    cout << "polygon:" << endl;
    std::copy(polygon.begin(), polygon.end(),
	std::ostream_iterator<Point>(cout, "\n"));
    cout << endl;

    iterator left =
	CGAL::left_vertex_2(polygon.begin(), polygon.end());
    iterator right =
	CGAL::right_vertex_2(polygon.begin(), polygon.end());
    iterator top =
	CGAL::top_vertex_2(polygon.begin(), polygon.end());
    iterator bottom =
	CGAL::bottom_vertex_2(polygon.begin(), polygon.end());
    bool simple =
	CGAL::is_simple_2(polygon.begin(), polygon.end());
    bool convex =
	CGAL::is_convex_2(polygon.begin(), polygon.end());
    CGAL::Bbox_2 bbox =
	CGAL::bbox_2(polygon.begin(), polygon.end(), R());
    typename R::FT area = 0, area2;
    CGAL::area_2(polygon.begin(), polygon.end(), area, R());
    area2 = CGAL::polygon_area_2(polygon.begin(), polygon.end(), R());
    CGAL::Bounded_side bside =
	 CGAL::bounded_side_2(polygon.begin(), polygon.end(), point);
    CGAL::Oriented_side oside =
	CGAL::oriented_side_2(polygon.begin(), polygon.end(), point);
    CGAL::Orientation orientation =
	CGAL::orientation_2(polygon.begin(), polygon.end());

    cout << "left   = " << *left << endl;
    cout << "right  = " << *right << endl;
    cout << "top    = " << *top << endl;
    cout << "bottom = " << *bottom << endl;
    cout << "the polygon is ";
    if (!simple) cout << "not ";
    cout << "simple" << endl;
    cout << "the polygon is ";
    if (!convex) cout << "not ";
    cout << "convex" << endl;
    cout << "the bounding box is " << bbox << endl;
    cout << "the area is " << area << endl;
    if (area != area2)
        cout << "but according to polygon_area_2 it is " << area2 << endl;

    switch (bside) {
    case CGAL::ON_BOUNDED_SIDE:
	cout << "the point is on bounded side" << endl;
	break;
    case CGAL::ON_BOUNDARY:
	cout << "the point is on the boundary" << endl;
	break;
    case CGAL::ON_UNBOUNDED_SIDE:
	cout << "the point is on the unbounded side" << endl;
	break;
    }

    switch (oside) {
    case CGAL::ON_NEGATIVE_SIDE:
	cout << "the point is on the negative side" << endl;
	break;
    case CGAL::ON_ORIENTED_BOUNDARY:
	cout << "the point is on the oriented boundary" << endl;
	break;
    case CGAL::ON_POSITIVE_SIDE:
	cout << "the point is on the positive side" << endl;
	break;
    }

    switch(orientation) {
    case CGAL::CLOCKWISE:
	cout << "the orientation is clockwise" << endl;
	break;
    case CGAL::COUNTERCLOCKWISE:
	cout << "the orientation is counter clockwise" << endl;
	break;
    case CGAL::COLLINEAR:
	cout << "the orientation is collinear" << endl;
	break;
    }
}

//-----------------------------------------------------------------------//
//                          main
//-----------------------------------------------------------------------//

int main()
{
  CGAL::set_pretty_mode(cout);

  cout << endl;
  cout << "--------------------------------------------------------" << endl;
  cout << "-   testing polygon algorithms with cartesian/double   -" << endl;
  cout << "--------------------------------------------------------" << endl;
  cout << endl;
  typedef CGAL::Cartesian<double> R1;

  typedef CGAL::Point_2<R1> Point1;
  test_polygon(R1(), Point1(), "data/polygon_cartesian.dat");

  cout << endl;
  cout << "--------------------------------------------------------" << endl;
  cout << "-   testing polygon algorithms with homogeneous/double -" << endl;
  cout << "--------------------------------------------------------" << endl;
  cout << endl;
  typedef CGAL::Homogeneous<double> R2;
  typedef CGAL::Point_2<R2> Point2;
  test_polygon(R2(), Point2(), "data/polygon_homogeneous.dat");

  return 0;
}

