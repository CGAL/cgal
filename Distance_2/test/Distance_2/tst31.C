#include "numrep1.h"
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/squared_distance_2_2.h>


#include "numrep2.h"

typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Line_2< TestR > line_t;
typedef CGAL::Triangle_2< TestR > triangle_t;


int main()
{
    randomint ri;
    int x1, x2, x3, y1, y2, y3, w1, w2, w3;
    TestR::FT d;
    std::cin >> x1 >> y1 >> x2 >> y2;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    line_t line(tp1, tp2);
    std::cin >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    w3 = ri.next();
    point_t tp3(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp4(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    point_t tp5(to_nt(w3*x3), to_nt(w3*y3), to_nt(w3));
    triangle_t tr(tp3, tp4, tp5);
    d = CGAL::squared_distance(line, tr);
    std::cout << CGAL::to_double(d) << '\n';
    return 0;
}
