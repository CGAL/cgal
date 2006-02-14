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

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>


typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Segment_2< TestR > segment_t;

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

void print(const segment_t &iseg)
{
	print(iseg.start());
	std::cout<<' ';
	print(iseg.end());
}


void treat_intersection(const segment_t &seg1, const segment_t &seg2)
{
    point_t pt1;
    segment_t iseg;
    bool is;
    is = CGAL::do_intersect(seg1,seg2);
    CGAL::Object result = CGAL::intersection(seg1, seg2);
    if (!CGAL::assign(iseg, result) && !CGAL::assign(pt1, result)) {
	std::cout << "No intersection.\n";
        assert(!is);
    }
    if (CGAL::assign(pt1, result)) {
	std::cout << "Point intersection.\n";
	print(pt1);
	std::cout<<'\n';
        assert(is);
    }
    if (CGAL::assign(iseg, result)) {
	std::cout << "Segment intersection.\n";
	print(iseg);
	std::cout<<'\n';
        assert(is);
    }
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
    point_t p1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1)),
	    p2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    segment_t seg1(p1, p2);
    std::cin >> x1 >> y1 >> x2 >> y2;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    p1 = point_t(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    p2 = point_t(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    segment_t seg2(p1, p2);
    treat_intersection(seg1, seg2);
    return 0;
}
