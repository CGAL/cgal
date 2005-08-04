#include "numrep1.h"
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/squared_distance_2.h>

#include "numrep2.h"

typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Line_2< TestR > line_t;



int main()
{
    randomint ri;
    int x1, x2, y1, y2, w1, w2;
    TestR::FT d;
    std::cin >> x1 >> y1 >> x2 >> y2;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    line_t line1(tp1, tp2);
    std::cin >> x1 >> y1 >> x2 >> y2;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp3(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp4(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    line_t line2(tp3, tp4);
    d = CGAL::squared_distance(line1, line2);
    std::cout << CGAL::to_double(d) << '\n';
    return 0;
}
