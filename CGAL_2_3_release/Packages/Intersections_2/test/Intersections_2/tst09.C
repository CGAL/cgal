#include "numrep1.h"
/*
 * Segment Ray intersection.
 */

#include <CGAL/Object.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Segment_2_Ray_2_intersection.h>




#include "numrep2.h"

typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Ray_2< TestR > ray_t;
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


void treat_intersection(const segment_t &seg, const ray_t &ray)
{
    point_t ipt;
    segment_t iseg;
    CGAL::Object result = CGAL::intersection(seg, ray);
    if (CGAL::assign(ipt, result)) {
	std::cout << "Point intersection.\n";
	print(ipt);
	std::cout<<'\n';
    }
    if (CGAL::assign(iseg, result)) {
	std::cout << "Segment intersection.\n";
	print(iseg.start());
	std::cout << ' ';
	print(iseg.end());
	std::cout<<'\n';
    }
    if (!CGAL::assign(iseg, result) && !CGAL::assign(ipt, result)) {
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
    w1 = ri.next();
    w2 = ri.next();
    point_t tp1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    segment_t seg(tp1, tp2);
    std::cin >> x1 >> y1 >> x2 >> y2;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp3(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp4(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    ray_t ray(tp3, tp4);
    treat_intersection(seg, ray);
    return 0;
}
