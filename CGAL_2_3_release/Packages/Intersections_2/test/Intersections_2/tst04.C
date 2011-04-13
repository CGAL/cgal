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
    std::cout.setf(std::ios::showpos, std::ios::showpos);
    std::cout << xd <<' '<< yd;
    std::cout.unsetf(std::ios::showpos);
}


void treat_intersection(const line_t &line1, const line_t &line2)
{
    
    point_t pt1;
    line_t l;
    CGAL::Object result = CGAL::intersection(line1, line2);
    if (CGAL::assign(pt1, result)) {
	std::cout << "Point intersection.\n";
	print(pt1);
	std::cout<<'\n';
    }
    if (CGAL::assign(l, result)) {
	std::cout << "Line intersection.\n";
    }
    if (!CGAL::assign(l, result) && !CGAL::assign(pt1, result))
	std::cout << "No intersection.\n";
}

int main()
{
    randomint ri;
    int x1, x2, y1, y2, w1, w2;
    std::cin >> x1 >> y1 >> x2 >> y2;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t p1, p2;
    p1 = point_t(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    p2 = point_t(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    line_t line1(p1, p2);
    std::cin >> x1 >> y1 >> x2 >> y2;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    p1 = point_t(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    p2 = point_t(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    line_t line2(p1, p2);
    treat_intersection(line1, line2);
    return 0;
}
