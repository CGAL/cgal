#include <CGAL/Cartesian.h>
#include "numrep1.h"

#include <iostream.h>
#include <CGAL/Object.h>
#include <CGAL/Triangle_2.h>
#include <vector.h>
#include <CGAL/Triangle_2_Triangle_2_intersection.h>

#include "numrep2.h"

typedef CGAL_Point_2< TestR > point_t;
typedef CGAL_Segment_2< TestR > segment_t;
typedef CGAL_Triangle_2< TestR > trian_t;
typedef vector<point_t>  container_t;

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

/*
template <class R>
void print(const CGAL_Point_2<R> &pt)
{
    cout << CGAL_to_double(pt.x()) <<' '<<CGAL_to_double(pt.y());
}
*/

void one_pair(trian_t const & trian1, trian_t const & trian2)
{
    segment_t seg;
    trian_t trian3;
    point_t point;
    container_t pgn;

    if (CGAL_do_intersect(trian1, trian2))
	;

    CGAL_Object result = CGAL_intersection(trian1, trian2);
    if (CGAL_assign(point, result)) {
	cout << "Point intersection.\n";
    }
    if (CGAL_assign(seg, result)) {
	cout << "Segment intersection.\n";
    }
    if (CGAL_assign(trian3, result)) {
	cout << "Triangle intersection.\n";
    }
    if (CGAL_assign(pgn, result)) {
	cout << "Polygon intersection.\n";
	cout << pgn.size()<<'\n';
	for (int i=0; i<pgn.size(); i++) {
	    print(pgn[i]);
	    cout << '\n';
	}
    }
    if (!CGAL_assign(seg, result) && !CGAL_assign(trian3, result)
	    && !CGAL_assign(point, result) && !CGAL_assign(pgn, result)) {
	cout << "No intersection.\n";
    }
}

int main()
{
    randomint ri;
    int x1, x2, x3, y1, y2, y3, w1, w2, w3;
    cin >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    w3 = ri.next();
    point_t tp1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    point_t tp3(to_nt(w3*x3), to_nt(w3*y3), to_nt(w3));
    trian_t trian1(tp1, tp2, tp3);
    cin >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;
    if (!cin)
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

/*
    CGAL_Triangle_2_Triangle_2_pair<TestR> pair(&trian1, &trian2);
    switch (pair.intersection_type()) {
    case CGAL_Triangle_2_Triangle_2_pair<TestR>::POLYGON:
	cout<<"Polygon intersection.\n";
	pair.intersection();
	break;
    case CGAL_Triangle_2_Triangle_2_pair<TestR>::TRIANGLE:
	cout<<"Triangle intersection.\n";
	pair.intersection(trian3);
	break;
    case CGAL_Triangle_2_Triangle_2_pair<TestR>::SEGMENT:
	cout<<"Segment intersection.\n";
	pair.intersection(seg);
	break;
    case CGAL_Triangle_2_Triangle_2_pair<TestR>::POINT:
	cout<<"Point intersection.\n";
	pair.intersection(point);
	break;
    case CGAL_Triangle_2_Triangle_2_pair<TestR>::NO:
	cout<<"No intersection.\n";
	break;
    }
*/
