#include <CGAL/Box_intersection_d.h>

#include <iostream>
#include <vector>
#include <CGAL/Random.h>

typedef CGAL::Box_intersection_d::Box_d< double, 3 > Box;

void fill_boxes( unsigned int n, std::vector<Box>& boxes ) {
    double min[3], max[3];
    for( unsigned int i = 0; i < n; ++i ) {
        for( unsigned int d = 0; d < 3; ++d ) {
            min[d] = 10 * CGAL::default_random.get_double();
            max[d] = min[d] + 1 + CGAL::default_random.get_double();
        }
        boxes.push_back( Box( min, max) );
    }
}

void callback( const Box& a, const Box& b ) {
    std::cout << "intersection between box "
              << a.get_id() << " and " << b.get_id() << std::endl;
};

int main() {
    std::vector<Box> a;
    fill_boxes( 100, a );
    CGAL::box_intersection_d( a.begin(), a.end(), callback );
}

