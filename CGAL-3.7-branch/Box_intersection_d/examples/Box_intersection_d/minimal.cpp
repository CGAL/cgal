#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_2.h>
#include <iostream>

typedef CGAL::Box_intersection_d::Box_d<double,2> Box;
typedef CGAL::Bbox_2                              Bbox;
                                                     // 9 boxes of a grid
Box boxes[9] = { Bbox( 0,0,1,1), Bbox( 1,0,2,1), Bbox( 2,0,3,1), // low
                 Bbox( 0,1,1,2), Bbox( 1,1,2,2), Bbox( 2,1,3,2), // middle
                 Bbox( 0,2,1,3), Bbox( 1,2,2,3), Bbox( 2,2,3,3)};// upper
// 2 selected boxes as query; center and upper right
Box query[2] = { Bbox( 1,1,2,2), Bbox( 2,2,3,3)};

void callback( const Box& a, const Box& b ) {
    std::cout << "box " << a.id() << " intersects box " << b.id() << std::endl;
}
int main() {
    CGAL::box_intersection_d( boxes, boxes+9, query, query+2, callback);
    return 0;
}
