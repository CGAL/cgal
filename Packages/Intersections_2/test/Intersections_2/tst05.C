/*
 * Intersection of two segments.
 */
#include "numrep1.h"

#include <CGAL/Object.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>
#include <CGAL/number_utils.h>

#include "numrep2.h"


typedef CGAL_Point_2< TestR > point_t;
typedef CGAL_Segment_2< TestR > segment_t;

void print(const point_t &pt)
{
    double xd = CGAL_to_double(pt.x());
    double yd = CGAL_to_double(pt.y());
    // force 0 to be positive zero.
    if (xd == 0.0)
	xd = 0.0;
    if (yd == 0.0)
	yd = 0.0;
    cout.setf(ios::showpos, ios::showpos);
    cout << xd <<' '<< yd;
    cout.setf(0, ios::showpos);
}

void print(const segment_t &iseg)
{
	print(iseg.start());
	cout<<' ';
	print(iseg.end());
}


void treat_intersection(const segment_t &seg1, const segment_t &seg2)
{
    point_t pt1;
    segment_t iseg;
    CGAL_Object result = CGAL_intersection(seg1, seg2);
    if (!CGAL_assign(iseg, result) && !CGAL_assign(pt1, result))
	cout << "No intersection.\n";
    if (CGAL_assign(pt1, result)) {
	cout << "Point intersection.\n";
	print(pt1);
	cout<<'\n';
    }
    if (CGAL_assign(iseg, result)) {
	cout << "Segment intersection.\n";
	print(iseg);
	cout<<'\n';
    }
}

int main()
{
    randomint ri;
    int x1, x2, y1, y2, w1, w2;
    cin >> x1 >> y1 >> x2 >> y2;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t p1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1)),
	    p2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    segment_t seg1(p1, p2);
    cin >> x1 >> y1 >> x2 >> y2;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    p1 = point_t(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    p2 = point_t(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    segment_t seg2(p1, p2);
    treat_intersection(seg1, seg2);
    return 0;
}
