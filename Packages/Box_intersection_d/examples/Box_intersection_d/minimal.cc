#include <CGAL/Box_intersection_d.h>

#include <iostream>
#include <cstdlib>
#include <vector>

typedef CGAL::Box_intersection_d::Box_d< double, 3 > Box;

void fill_boxes( unsigned int n, std::vector<Box>& boxes ) {
    double lo[3], hi[3];
    for( unsigned int i = 0; i < n; ++i ) {
        for( unsigned int d = 0; d < 3; ++d ) {
            lo[d] = 10 * drand48();
            hi[d] = lo[d] + 1 + drand48();
        }
        boxes.push_back( Box( lo, hi) );
    }
}

void callback( const Box& a, const Box& b ) {
    std::cout << "intersection between box "
              << a.get_id() << " and " << b.get_id() << std::endl;
};

int main() {
    std::vector<Box> a, b;
    fill_boxes( 100, a );
    fill_boxes( 100, b );
    CGAL::box_intersection_d( a.begin(), a.end(),
                              b.begin(), b.end(), callback );
}

