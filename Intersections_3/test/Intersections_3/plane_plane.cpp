#include "numrep1.h"
#include <CGAL/Plane_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/intersection_3_1.h>
#include "numrep2.h"

using std::cout;
using std::cin;

typedef CGAL::Plane_3<TestR> Plane;
typedef CGAL::Line_3<TestR> Line;
typedef CGAL::Point_3<TestR> Point;
typedef CGAL::Direction_3<TestR> Direction;
typedef CGAL::Object  Object;


bool read_data(Plane &pl1, Plane &pl2)
{
    randomint ri;
    int a, b, c, d, w;
    cin >> a >> b >> c >> d;
    if (!cin)
	return false;
    w = ri.next();
    pl1 = Plane(a*w, b*w, c*w, d*w);
    cin >> a >> b >> c >> d;
    if (!cin)
	return false;
    w = ri.next();
    pl2 = Plane(a*w, b*w, c*w, d*w);
    return true;
}

void write_point(const Point & pt)
{
    double xd = CGAL::to_double(pt.x());
    double yd = CGAL::to_double(pt.y());
    double zd = CGAL::to_double(pt.z());
    // force 0 to be positive zero.
    if (xd == 0.0)
	xd = 0.0;
    if (yd == 0.0)
	yd = 0.0;
    if (zd == 0.0)
	zd = 0.0;
    cout << xd <<' '<< yd <<' '<< zd << '\n';
}

// Testing intersections of 3 planes.
bool plane_plane_plane()
{
  Plane pl1(1,0,0,0);
  Plane pl2(0,1,0,0);
  Plane pl3(0,0,1,0);

  // Generic intersection.
  Object o = CGAL::intersection(pl1, pl2, pl3);
  Point p;
  if (!assign(p, o)) {
    std::cerr << "Unexpected intersection result" << std::endl;
    return false;
  }

  if (p != Point(0,0,0)) {
    std::cerr << "Unexpected intersection result" << std::endl;
    return false;
  }

  // Empty intersection.
  Plane pl4(1,0,0,1); // pl4 is // to pl1.

  Object o2 = CGAL::intersection(pl1, pl2, pl4);
  if (!o2.is_empty()) {
    std::cerr << "Unexpected intersection result" << std::endl;
    return false;
  }

  Object o3 = CGAL::intersection(pl1, pl4, pl2);
  if (!o3.is_empty()) {
    std::cerr << "Unexpected intersection result" << std::endl;
    return false;
  }

  // Intersection in a line.
  Plane pl5(1,1,0,0); // pl1, pl2, pl5 intersect in the line l.
  Line l;

  Object o4 = CGAL::intersection(pl1, pl2, pl5);
  if (!assign(l, o4)) {
    std::cerr << "Unexpected intersection result" << std::endl;
    return false;
  }

  if (l != Line(Point(0,0,0), Point(0,0,1))) {
    std::cerr << "Unexpected intersection result" << std::endl;
    return false;
  }

  // Intersection in a plane.
  Object o5 = CGAL::intersection(pl1, pl1, pl1);
  Plane pl;
  if (!assign(pl, o5)) {
    std::cerr << "Unexpected intersection result" << std::endl;
    return false;
  }

  if (pl != pl1) {
    std::cerr << "Unexpected intersection result" << std::endl;
    return false;
  }

  return true;
}

int main()
{
    if (!plane_plane_plane())
      return 1;

    Plane pl1, pl2;
    Line l;
    CGAL::Object result;
    
    if (!read_data(pl1, pl2))
	return 1;
    result = CGAL::intersection(pl1, pl2);
    if (result.is_empty()) {
	if (CGAL::do_intersect(pl1, pl2)) {
	    cout << "do_intersect is inconsistent with intersection result.\n";
	    return 1;
	}
	cout << "No intersection.\n";
	return 0;
    }
    if (!CGAL::do_intersect(pl1, pl2)) {
	cout << "do_intersect is inconsistent with intersection result.\n";
	return 1;
    }
    if (CGAL::assign(l, result)) {
	cout << "Line intersection.\n";
	write_point(l.point(0));
	write_point(l.point(1));
	return 0;
    } else {
	if (CGAL::assign(pl1, result)) {
	    cout << "Plane intersection.\n";
	    return 0;
	}
    }
    cout << "Unknown result.\n";
    return 1;
}
