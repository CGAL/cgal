#include "numrep1.h"
#include <CGAL/Point_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/squared_distance_2.h>


#include "numrep2.h"

typedef CGAL_Point_2< TestR > point_t;
typedef CGAL_Ray_2< TestR > ray_t;


int main()
{
    randomint ri;
    int x1, x2, y1, y2, w1, w2;
    TestR::FT d;
    cin >> x1 >> y1;
    if (!cin)
	return 1;
    w1 = ri.next();
    point_t pt(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    cin >> x1 >> y1 >> x2 >> y2;
    if (!cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    ray_t ray(tp1, tp2);
    d = CGAL_squared_distance(pt, ray);
    cout << CGAL_to_double(d) << '\n';
    return 0;
}
