/*
 * 2D Iso_rectangle Ray intersection.
 */

#include "numrep1.h"
#include <CGAL/Object.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Ray_2_Iso_rectangle_2_intersection.h>
#include "numrep2.h"
#include <iostream>

using std::cout;
using std::cin;
using std::ios;

typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Ray_2< TestR > ray_t;
typedef CGAL::Segment_2< TestR > segment_t;
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
    cout.setf(ios::showpos, ios::showpos);
    cout << xd <<' '<< yd;
    cout.setf(0, ios::showpos);
}



void one_pair(ray_t const & ray, rect_t const & rect)
{
    segment_t seg;
    point_t point;
/*
    CGAL::Ray_2_Iso_rectangle_2_pair<TestR> pair(&ray, &rect);
    switch (pair.intersection_type()) {
    case CGAL::Ray_2_Iso_rectangle_2_pair<TestR>::SEGMENT:
	cout<<"Segment intersection.\n";
	pair.intersection(seg);
	print(seg.start());
	cout<<' ';
	print(seg.end());
	cout<<'\n';
	break;
    case CGAL::Ray_2_Iso_rectangle_2_pair<TestR>::POINT:
	cout<<"Point intersection.\n";
	pair.intersection(point);
	print(point);
	cout<<'\n';
	break;
    case CGAL::Ray_2_Iso_rectangle_2_pair<TestR>::NO:
	cout<<"No intersection.\n";
	break;
    }
*/
    CGAL::Object result = CGAL::intersection(ray, rect);
    if (CGAL::assign(point, result)) {
	cout << "Point intersection.\n";
    }
    if (CGAL::assign(seg, result)) {
	cout << "Segment intersection.\n";
    }
    if (!CGAL::assign(point, result) && !CGAL::assign(seg, result)) {
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
    point_t p1, p2;
    w1 = ri.next();
    w2 = ri.next();
    p1 = point_t(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    p2 = point_t(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    ray_t ray(p1, p2);
    cin >> x1 >> y1 >> x2 >> y2;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    p1 = point_t(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    p2 = point_t(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    rect_t rect(p1, p2);
    one_pair(ray, rect);
    return 0;
}
