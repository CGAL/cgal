#include "numrep1.h"
#include <CGAL/Object.h>
#include <CGAL/number_utils.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Line_2_Ray_2_intersection.h>

#include "numrep2.h"

typedef CGAL_Point_2< TestR > point_t;
typedef CGAL_Line_2< TestR > line_t;
typedef CGAL_Ray_2< TestR > ray_t;
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


void treat_intersection(const line_t &line, const ray_t &ray)
{
    point_t ipt;
    ray_t iray;
/*
    CGAL_Line_2_Ray_2_pair<TestR> pair(&line, &ray);
    switch (pair.intersection_type()) {
    case CGAL_Line_2_Ray_2_pair<TestR>::NO:
	cout << "No intersection.\n";
	break;
    case CGAL_Line_2_Ray_2_pair<TestR>::RAY:
	cout << "Ray intersection.\n";
	pair.intersection(iray);
	print(iray.start());
	cout << ' ';
	print(CGAL_ORIGIN + iray.direction().vector());
	cout<<'\n';
	break;
    case CGAL_Line_2_Ray_2_pair<TestR>::POINT:
	cout << "Point intersection.\n";
	pair.intersection(ipt);
	print(ipt);
	cout<<'\n';
	break;
    default:
	cout << "Unexpected result.\n";
    }
*/

    CGAL_Object result = CGAL_intersection(line, ray);
    if (CGAL_assign(ipt, result)) {
	cout << "Point intersection.\n";
	print(ipt);
	cout<<'\n';
    }
    if (CGAL_assign(iray, result)) {
	cout << "Ray intersection.\n";
	print(iray.start());
	cout << ' ';
	print(CGAL_ORIGIN + iray.direction().vector());
	cout<<'\n';
    }
    if (!CGAL_assign(iray, result) && !CGAL_assign(ipt, result)) {
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
    line_t line(tp1, tp2);
    cin >> x1 >> y1 >> x2 >> y2;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp3(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp4(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    ray_t ray(tp3, tp4);
    treat_intersection(line, ray);
    return 0;
}
