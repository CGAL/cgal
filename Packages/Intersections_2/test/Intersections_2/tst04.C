/*
 * Line Line intersection.
 */

#include "numrep1.h"
#include <CGAL/Object.h>
#include <CGAL/Line_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Line_2_Line_2_intersection.h>


#include "numrep2.h"

typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Line_2< TestR > line_t;

void print(const point_t &pt)
{
    double xd = CGAL::to_double(pt.x());
    double yd = CGAL::to_double(pt.y());
    // force 0 to be positive zero.
    if (xd == 0.0)
	xd = 0.0;
    if (yd == 0.0)
	yd = 0.0;
    cout.setf(ios::showpos, ios::showpos);
    cout << xd <<' '<< yd;
    cout.setf(0, ios::showpos);
}


void treat_intersection(const line_t &line1, const line_t &line2)
{
    
    point_t pt1;
    line_t l;
/*
    typedef CGAL::Line_2_Line_2_pair<TestR> is_t;
    is_t linepair(&line1, &line2);
    switch ( linepair.intersection_type()) {
    case is_t::NO:
	cout << "No intersection.\n";
        break;
    case is_t::POINT:
	cout << "Point intersection.\n";
        linepair.intersection(pt1);
	print(pt1);
	cout<<'\n';
        break;
    case is_t::LINE:
	cout << "Line intersection.\n";
        linepair.intersection(l);
        break;
    default:
	cout << "Unexpected result.\n";
    }
*/
    CGAL::Object result = CGAL::intersection(line1, line2);
    if (CGAL::assign(pt1, result)) {
	cout << "Point intersection.\n";
	print(pt1);
	cout<<'\n';
    }
    if (CGAL::assign(l, result)) {
	cout << "Line intersection.\n";
    }
    if (!CGAL::assign(l, result) && !CGAL::assign(pt1, result))
	cout << "No intersection.\n";
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
    point_t p1, p2;
    p1 = point_t(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    p2 = point_t(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    line_t line1(p1, p2);
    cin >> x1 >> y1 >> x2 >> y2;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    p1 = point_t(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    p2 = point_t(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    line_t line2(p1, p2);
    treat_intersection(line1, line2);
    return 0;
}
