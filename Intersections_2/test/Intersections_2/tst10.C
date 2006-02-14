#include "numrep1.h"
#include <CGAL/Object.h>
#include <CGAL/number_utils.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Line_2_Ray_2_intersection.h>

#include "numrep2.h"

typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Line_2< TestR > line_t;
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


void treat_intersection(const line_t &line, const ray_t &ray)
{
    point_t ipt;
    ray_t iray;

    CGAL::Object result = CGAL::intersection(line, ray);
    if (CGAL::assign(ipt, result)) {
	std::cout << "Point intersection.\n";
	print(ipt);
	std::cout<<'\n';
    }
    if (CGAL::assign(iray, result)) {
	std::cout << "Ray intersection.\n";
	print(iray.start());
	std::cout << ' ';
	print(CGAL::ORIGIN + iray.direction().vector());
	std::cout<<'\n';
    }
    if (!CGAL::assign(iray, result) && !CGAL::assign(ipt, result)) {
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
    line_t line(tp1, tp2);
    std::cin >> x1 >> y1 >> x2 >> y2;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp3(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp4(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    ray_t ray(tp3, tp4);
    treat_intersection(line, ray);
    return 0;
}
