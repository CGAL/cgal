#include <CGAL/Box_intersection_d.h>

#include <iostream>
#include <cstdlib>
#include <vector>

using namespace std;
using namespace CGAL;

struct Primitive {
    float f;
    Primitive() : f( drand48() ) {}
    bool intersect( Primitive *other )
    { return f * other->f > 0.6f; }
};

struct Box : public UniqueNumbers
{
    Primitive *primitive;
    float lo[3], size;
    Box( Primitive *p ) : primitive( p ) {
        for( unsigned int d = 0; d < 3; ++d )
            lo[d] = 10.0 * drand48();
        size = 1 + 3.0 * drand48();
    }
};

struct BoxAdapter {
    typedef ::Box Box;
    typedef float NumberType;

    static NumberType get_lo( const Box& b, unsigned int dim )
    { return b.lo[ dim ]; }

    static NumberType get_hi( const Box& b, unsigned int dim )
    { return b.lo[ dim ] + b.size; }

    static unsigned int get_num( const Box& b )
    { return b.num();     }

    static unsigned int get_dim() { return 3; }
};

typedef vector< Box > BoxContainer;
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
                  boxes2.begin(), boxes2.end(), callback, Traits(), 2 );
}

