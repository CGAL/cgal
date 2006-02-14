#include "numrep1.h"
#include <assert.h>
/*
 * Squared distance between two segments.
 */
#include <CGAL/Point_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/squared_distance_3.h>

#include "numrep2.h"

typedef CGAL::Point_3< TestR > point_t;
typedef CGAL::Segment_3< TestR > segment_t;
typedef CGAL::Ray_3< TestR > ray_t;

using std::cin;
using std::cout;

int main()
{
    randomint ri;
    int x1, x2, y1, y2, z1, z2, w1, w2;
    TestR::FT d;
    cin >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1*z1), to_nt(w1));
    point_t tp2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2*z2), to_nt(w2));
    segment_t seg(tp1, tp2);
    cin >> x1 >> y1 >> z1 >> x2 >> y2 >> z2;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp3(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1*z1), to_nt(w1));
    point_t tp4(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2*z2), to_nt(w2));
    ray_t ray(tp3, tp4);
    d = CGAL::squared_distance(seg, ray);
    cout << CGAL::to_double(d) << '\n';
    return 0;
}
