#include <CGAL/Box_intersection_d/box_traits.h>
#include <CGAL/Box_intersection_d/segment_tree.h>

#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;

typedef Default_Bbox_d< double, 3 >      Box;
typedef Default_Bbox_d_Adapter< Box >    BoxAdapter;
typedef Default_Box_Traits< BoxAdapter > Traits;
typedef vector< Box >                    BoxContainer;

void fill_boxes( unsigned int n, BoxContainer& boxes ) {
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
    cout << "intersection between box "
         << a.num() << " and " << b.num() << endl;
};

int main() {
    BoxContainer boxes1, boxes2;
    fill_boxes( 100, boxes1 );
    fill_boxes( 100, boxes2 );
    segment_tree( boxes1.begin(), boxes1.end(),
                  boxes2.begin(), boxes2.end(), callback, Traits(), 2 );
}

