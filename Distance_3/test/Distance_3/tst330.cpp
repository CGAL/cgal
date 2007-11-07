#include "numrep1.h"
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/squared_distance_3_2.h>


#include "numrep2.h"

typedef CGAL::Point_3< TestR > point_t;
typedef CGAL::Plane_3< TestR > plane_t;

using std::cin;
using std::cout;

int main()
{
    randomint ri;
    int x, y, z, w, w1;
    TestR::FT d;
    cin >> x >> y >> z;
    if (!cin)
	return 1;
    w1 = ri.next();
    point_t pt(to_nt(w1*x), to_nt(w1*y), to_nt(w1*z), to_nt(w1));
    cin >> x >> y >> z >> w;
    if (!cin)
	return 1;
    plane_t plane(to_nt(x), to_nt(y), to_nt(z), to_nt(w));
    d = CGAL::squared_distance(pt, plane);
    cout << CGAL::to_double(d) << '\n';
    return 0;
}
