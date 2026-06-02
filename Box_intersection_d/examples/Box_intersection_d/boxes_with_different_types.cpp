#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Box_intersection_d/Box_d.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel           Kernel;
typedef Kernel::Point_2                                               Point_2;
typedef Kernel::Segment_2                                             Segment_2;
typedef CGAL::Box_intersection_d::Box_with_info_d<double,2,Point_2>   Point_box;
typedef CGAL::Box_intersection_d::Box_with_info_d<double,2,Segment_2> Segment_box;
                                                     // 4 points of a grid
Point_box pts[4] = { Point_box( Point_2(0, 0).bbox(), Point_2(0, 0)),
                     Point_box( Point_2(0, 1).bbox(), Point_2(0, 1)),
                     Point_box( Point_2(1, 0).bbox(), Point_2(1, 0)),
                     Point_box( Point_2(1, 1).bbox(), Point_2(1, 1))};// upper
// 2 segment as query
Segment_2 seg1(Point_2(0.2,-0.2), Point_2(1.2, 0.8));
Segment_2 seg2(Point_2(-0.2,1.3), Point_2(1.2,-0.1));
Segment_box query[2] = { Segment_box(seg1.bbox(), seg1),
                         Segment_box(seg2.bbox(), seg2)};

void callback( const Point_box& a, const Segment_box& b ) {
    std::cout << "box of point " << a.info() << " intersects box of segment " << b.info() << std::endl;
}
int main() {
    CGAL::box_intersection_d( pts, pts+4, query, query+2, callback);
    return 0;
}
