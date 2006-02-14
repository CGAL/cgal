#include "numrep1.h"
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/squared_distance_2_1.h>


#include "numrep2.h"

typedef CGAL::Point_2< TestR > point_t;
typedef CGAL::Line_2< TestR > line_t;


int main()
{
    randomint ri;
    int x1, x2, y1, y2, w1, w2;
    TestR::FT d,d2;
    std::cin >> x1 >> y1;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    point_t pt(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    std::cin >> x1 >> y1 >> x2 >> y2;
    if (!std::cin)
	return 1;
    w1 = ri.next();
    w2 = ri.next();
    point_t tp1(to_nt(w1*x1), to_nt(w1*y1), to_nt(w1));
    point_t tp2(to_nt(w2*x2), to_nt(w2*y2), to_nt(w2));
    line_t line(tp1, tp2);
    d = CGAL::squared_distance(pt, line);
    std::cout << CGAL::to_double(d) << '\n';
    d2 = CGAL::Squared_distance_to_line<TestR>(line)(pt);
    if (d != d2)
        std::cout<<"Two computation methods give different results: "
                 << d << " and "<< d2 <<"\n";
    return (d2 == d) ? 0 : 1;
}
