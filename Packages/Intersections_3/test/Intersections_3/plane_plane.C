#include "numrep1.h"
#include <CGAL/Plane_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/intersection_3_1.h>
#include "numrep2.h"

typedef CGAL_Plane_3<TestR> Plane;
typedef CGAL_Line_3<TestR> Line;
typedef CGAL_Point_3<TestR> Point;
typedef CGAL_Direction_3<TestR> Direction;


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
    double xd = CGAL_to_double(pt.x());
    double yd = CGAL_to_double(pt.y());
    double zd = CGAL_to_double(pt.z());
    // force 0 to be positive zero.
    if (xd == 0.0)
	xd = 0.0;
    if (yd == 0.0)
	yd = 0.0;
    if (zd == 0.0)
	zd = 0.0;
    cout << xd <<' '<< yd <<' '<< zd << '\n';
}


int main()
{
    Plane pl1, pl2;
    Line l;
    CGAL_Object result;
    
    if (!read_data(pl1, pl2))
	return 1;
    result = CGAL_intersection(pl1, pl2);
    if (CGAL_assign(l, result)) {
	cout << "Line intersection.\n";
	write_point(l.point(0));
	write_point(l.point(1));
    } else {
	if (CGAL_assign(pl1, result)) {
	    cout << "Plane intersection.\n";
	} else {
	    cout << "No intersection.\n";
	}
    }
    return 0;
}
