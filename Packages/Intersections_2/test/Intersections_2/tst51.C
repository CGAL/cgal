#include "numrep1.h"
#include <CGAL/Object.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Line_2_Triangle_2_intersection.h>
#include <iostream.h>


#include "numrep2.h"

typedef CGAL_Point_2< TestR > point_t;
typedef CGAL_Line_2< TestR > line_t;
typedef CGAL_Segment_2< TestR > segment_t;
typedef CGAL_Triangle_2< TestR > trian_t;

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



void one_pair(line_t const & line, trian_t const & trian)
{
    segment_t seg;
    point_t point;
/*
    CGAL_Line_2_Triangle_2_pair<TestR> pair(&line, &trian);
    switch (pair.intersection_type()) {
    case CGAL_Line_2_Triangle_2_pair<TestR>::SEGMENT:
	cout<<"Segment intersection.\n";
	pair.intersection(seg);
	print(seg.start());
	cout<<' ';
	print(seg.end());
	cout<<'\n';
	break;
    case CGAL_Line_2_Triangle_2_pair<TestR>::POINT:
	cout<<"Point intersection.\n";
	pair.intersection(point);
	print(point);
	cout<<'\n';
	break;
    case CGAL_Line_2_Triangle_2_pair<TestR>::NO:
	cout<<"No intersection.\n";
	break;
    }
*/

    if (CGAL_do_intersect(line, trian))
	;

    CGAL_Object result = CGAL_intersection(line, trian);
    if (CGAL_assign(point, result)) {
	cout << "Point intersection.\n";
    }
    if (CGAL_assign(seg, result)) {
	cout << "Segment intersection.\n";
    }
    if (!CGAL_assign(seg, result) && !CGAL_assign(point, result)) {
	cout << "No intersection.\n";
    }
}

int main()
{
    randomint ri;
    int x1, x2, x3, y1, y2, y3, w1, w2, w3;
    cin >> x1 >> y1 >> x2 >> y2;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    line_t line(tp1, tp2);
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
    one_pair(line, trian);
    return 0;
}
