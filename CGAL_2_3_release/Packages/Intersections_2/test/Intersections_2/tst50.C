#include <CGAL/Cartesian.h>
#include "numrep1.h"

#include <iostream>
#include <CGAL/Object.h>
#include <CGAL/Triangle_2.h>
#include <vector>
#include <CGAL/Triangle_2_Triangle_2_intersection.h>

#include "numrep2.h"

typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Segment_2< TestR > segment_t;
typedef CGAL::Triangle_2< TestR > trian_t;
typedef std::vector<point_t>  container_t;

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

void one_pair(trian_t const & trian1, trian_t const & trian2)
{
    segment_t seg;
    trian_t trian3;
    point_t point;
    container_t pgn;

    if (CGAL::do_intersect(trian1, trian2))
	;

    CGAL::Object result = CGAL::intersection(trian1, trian2);
    if (CGAL::assign(point, result)) {
	std::cout << "Point intersection.\n";
    }
    if (CGAL::assign(seg, result)) {
	std::cout << "Segment intersection.\n";
    }
    if (CGAL::assign(trian3, result)) {
	std::cout << "Triangle intersection.\n";
    }
    if (CGAL::assign(pgn, result)) {
	std::cout << "Polygon intersection.\n";
	std::cout << pgn.size()<<'\n';
	for (container_t::size_type i=0; i<pgn.size(); i++) {
	    print(pgn[i]);
	    std::cout << '\n';
	}
    }
    if (!CGAL::assign(seg, result) && !CGAL::assign(trian3, result)
	    && !CGAL::assign(point, result) && !CGAL::assign(pgn, result)) {
	std::cout << "No intersection.\n";
    }
}

int main()
{
    randomint ri;
    int x1, x2, x3, y1, y2, y3, w1, w2, w3;
    std::cin >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    w3 = ri.next();
    point_t tp1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    point_t tp3(to_nt(w3*x3), to_nt(w3*y3), to_nt(w3));
    trian_t trian1(tp1, tp2, tp3);
    std::cin >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    w3 = ri.next();
    point_t tp4(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp5(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    point_t tp6(to_nt(w3*x3), to_nt(w3*y3), to_nt(w3));
    trian_t trian2(tp4, tp5, tp6);
    one_pair(trian1, trian2);
    return 0;
}
