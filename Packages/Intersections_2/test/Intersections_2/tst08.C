#include "numrep1.h"
#include <CGAL/Object.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Ray_2_Ray_2_intersection.h>

#include "numrep2.h"

using std::cout;
using std::cin;
using std::ios;

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
    cout.setf(ios::showpos, ios::showpos);
    cout << xd <<' '<< yd;
    cout.setf(0, ios::showpos);
}


void treat_intersection(const ray_t &ray1, const ray_t &ray2)
{
    point_t ipt;
    ray_t iray;
    segment_t iseg;
/*
    typedef CGAL::Ray_2_Ray_2_pair<TestR> is_t;
    is_t pair(&ray1, &ray2);
    switch (pair.intersection_type()) {
    case is_t::NO:
	cout << "No intersection.\n";
	break;
    case is_t::RAY:
	cout << "Ray intersection.\nstart: ";
	pair.intersection(iray);
	print(iray.start());
	cout<<" direction: ";
	print(CGAL::ORIGIN + iray.direction().vector());
	cout<<'\n';
	break;
    case is_t::SEGMENT:
	cout << "Segment intersection.\n";
	pair.intersection(iseg);
	print(iseg.start());
	cout << ' ';
	print(iseg.end());
	cout<<'\n';
	break;
    case is_t::POINT:
	cout << "Point intersection.\n";
	pair.intersection(ipt);
	print(ipt);
	cout<<'\n';
	break;
    default:
	cout << "Unexpected result.\n";
    }
*/

    CGAL::Object result = CGAL::intersection(ray1, ray2);
    if (CGAL::assign(ipt, result)) {
	cout << "Point intersection.\n";
	print(ipt);
	cout<<'\n';
    }
    if (CGAL::assign(iseg, result)) {
	cout << "Segment intersection.\n";
	print(iseg.start());
	cout << ' ';
	print(iseg.end());
	cout<<'\n';
    }
    if (CGAL::assign(iray, result)) {
	cout << "Ray intersection.\n";
	print(iray.start());
	cout<<" direction: ";
	print(CGAL::ORIGIN + iray.direction().vector());
	cout<<'\n';
    }
    if (!CGAL::assign(iseg, result) && !CGAL::assign(iray, result)
	    && !CGAL::assign(ipt, result)) {
	cout << "No intersection.\n";
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
    point_t tp1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    ray_t ray1(tp1, tp2);
    cin >> x1 >> y1 >> x2 >> y2;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp3(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp4(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    ray_t ray2(tp3, tp4);
    treat_intersection(ray1, ray2);
    return 0;
}
