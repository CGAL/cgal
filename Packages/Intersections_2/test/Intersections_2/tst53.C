/*
 * 2D Triangle Segment intersection.
 */

#include "numrep1.h"

#include <CGAL/Object.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Segment_2_Triangle_2_intersection.h>

#include "numrep2.h"
#include <iostream.h>

typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Segment_2< TestR > segment_t;
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
    cout.setf(ios::showpos, ios::showpos);
    cout << xd <<' '<< yd;
    cout.setf(0, ios::showpos);
}



void one_pair(segment_t const & seg, trian_t const & trian)
{
    segment_t iseg;
    point_t point;
/*
    CGAL::Segment_2_Triangle_2_pair<TestR> pair(&seg, &trian);
    switch (pair.intersection_type()) {
    case CGAL::Segment_2_Triangle_2_pair<TestR>::SEGMENT:
	cout<<"Segment intersection.\n";
	pair.intersection(iseg);
	print(iseg.start());
	cout<<' ';
	print(iseg.end());
	cout<<'\n';
	break;
    case CGAL::Segment_2_Triangle_2_pair<TestR>::POINT:
	cout<<"Point intersection.\n";
	pair.intersection(point);
	print(point);
	cout<<'\n';
	break;
    case CGAL::Segment_2_Triangle_2_pair<TestR>::NO:
	cout<<"No intersection.\n";
	break;
    }
*/

    if (CGAL::do_intersect(seg, trian))
	;

    CGAL::Object result = CGAL::intersection(seg, trian);
    if (CGAL::assign(point, result)) {
	cout << "Point intersection.\n";
    }
    if (CGAL::assign(iseg, result)) {
	cout << "Segment intersection.\n";
    }
    if (!CGAL::assign(iseg, result) && !CGAL::assign(point, result)) {
	cout << "No intersection.\n";
    }
}

main()
{
    int x1, x2, x3, y1, y2, y3, w1, w2, w3;
    randomint ri;
    cin >> x1 >> y1 >> x2 >> y2;
    if (!cin)
	return 1;
    
    w1 = ri.next();
    w2 = ri.next();
    point_t tp1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    segment_t seg(tp1, tp2);
    cin >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    w3 = ri.next();
    point_t tp3(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp4(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    point_t tp5(to_nt(w3*x3), to_nt(w3*y3), to_nt(w3));
    trian_t trian(tp3, tp4, tp5);
    one_pair(seg, trian);
    return 0;
}
