#include <CGAL/Cartesian.h>
#define CGAL_POLYGON_2_CONST_ITER
#include <CGAL/Polygon_2.h>

using CGAL::Polygon_2;
typedef CGAL::Cartesian<double> K;

int main()
{
    typedef K::Point_2 Point;
    Point a[] = {Point(0,0), Point(1,0), Point(1,1)};
    const int N = sizeof(a)/sizeof(a[0]);
    Polygon_2<K> pgn(a, a+N);
    *pgn.vertices_begin() = Point(0,1);
    return 0;
}

