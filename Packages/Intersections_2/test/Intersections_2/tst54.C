/*
 * 2D Triangle Point intersection.
 */

#include "numrep1.h"
#include <iostream>
#include <CGAL/Object.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Point_2_Triangle_2_intersection.h>

#include "numrep2.h"

typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Triangle_2< TestR > trian_t;

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



void one_pair(point_t const & pt, trian_t const & trian)
{
    point_t point;

    if (CGAL::do_intersect(pt, trian))
	;

    CGAL::Object result = CGAL::intersection(pt, trian);
    if (CGAL::assign(point, result)) {
	std::cout << "Point intersection.\n";
    }
    if (!CGAL::assign(point, result)) {
	std::cout << "No intersection.\n";
    }
}

int main()
{
    randomint ri;
    int x1, x2, x3, y1, y2, y3, w1, w2, w3;
    std::cin >> x1 >> y1;
    if (!std::cin)
	return 1;
    point_t p1, p2, p3;
    w1 = ri.next();
    point_t pt(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    std::cin >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    w3 = ri.next();
    p1 = point_t(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    p2 = point_t(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    p3 = point_t(to_nt(w3*x3), to_nt(w3*y3), to_nt(w3));
    trian_t trian(p1, p2, p3);
    one_pair(pt, trian);
    return 0;
}
