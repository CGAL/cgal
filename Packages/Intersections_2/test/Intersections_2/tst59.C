/*
 * 2D Iso_rectangle Segment intersection.
 */

#include "numrep1.h"
#include <CGAL/Object.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Iso_rectangle_2_Iso_rectangle_2_intersection.h>

#include "numrep2.h"
#include <iostream>

typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Iso_rectangle_2< TestR > rect_t;

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



void one_pair(rect_t const & irect1, rect_t const & irect2)
{
    rect_t irect;
    CGAL::Object result = CGAL::intersection(irect1, irect2);
    if (CGAL::assign(irect, result)) {
	std::cout << "Intersection.\n";
	print(irect.min());
	std::cout << ' ';
	print(irect.max());
	std::cout << '\n';
    }
    if (!CGAL::assign(irect, result)) {
	std::cout << "No intersection.\n";
    }
}

int main()
{
    randomint ri;
    int x1, x2, y1, y2, w1, w2;
    std::cin >> x1 >> y1 >> x2 >> y2;
    if (!std::cin)
	return 1;
    point_t p1, p2;
    w1 = ri.next();
    w2 = ri.next();
    p1 = point_t(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    p2 = point_t(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    rect_t irect1(p1, p2);
    std::cin >> x1 >> y1 >> x2 >> y2;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    p1 = point_t(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    p2 = point_t(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    rect_t irect2(p1, p2);
    one_pair(irect1, irect2);
    return 0;
}
