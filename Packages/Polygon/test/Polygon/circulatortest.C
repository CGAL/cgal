#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/circulator.h>
#include <CGAL/Polygon_2.h>
#include <list>

typedef CGAL::Cartesian<double> Rep;
typedef CGAL::Polygon_traits_2<Rep> Traits;
typedef CGAL::Point_2<Rep> Point;
typedef CGAL::Polygon_2<Traits, list<Point> > Polygon;
typedef CGAL::Circulator_from_container< std::vector<Point> >  Circulator;

int main()
{
    std::vector<Point> pts;
    Circulator c1(&pts);
    Polygon p1(c1);
    pts.push_back(Point(-3.1, 1.0));
    pts.push_back(Point(3.1, 1.0));
    pts.push_back(Point(1.1, 2.0));
    Circulator c2(&pts);
    Polygon p2(c2);
    return 0;
}
