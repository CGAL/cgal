#include <CGAL/Box_intersection_d.h>

#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;
using namespace CGAL;

struct Primitive {
    double f;
    Primitive() : f( drand48() ) {}
    bool intersect( Primitive *other )
    { return f * other->f > 0.6; }
};

struct Box : public Default_Bbox_d< double, 3 >
{
    Primitive *primitive;
    Box( Primitive *p ) : primitive( p ) {
        for( unsigned int d = 0; d < 3; ++d ) {
            _lo[d] = 10.0 * drand48();
            _hi[d] = _lo[d] + 1.0 + drand48();
        }
    }
};

typedef vector< Box >                  BoxContainer;
typedef Default_Bbox_d_Adapter< Box >  BoxAdapter;
typedef Default_Box_Traits< BoxAdapter, true > Traits;

void fill_boxes( unsigned int n, BoxContainer &boxes ) {
    for( unsigned int i = 0; i < n; ++i )
        boxes.push_back( Box( new Primitive() ) );
}

void callback( const Box &a, const Box &b ) {
    if( a.primitive->intersect( b.primitive ) )
        cout << "intersection between box "
            << a.num() << " and " << b.num() << endl;
};

int main() {
    BoxContainer boxes1, boxes2;
    fill_boxes( 100, boxes1 );
    fill_boxes( 100, boxes2 );
    segment_tree( boxes1.begin(), boxes1.end(),
                  boxes2.begin(), boxes2.end(), callback, Traits() );
}

