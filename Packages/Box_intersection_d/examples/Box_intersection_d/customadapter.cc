#include <CGAL/Box_intersection_d.h>

#include <iostream>
#include <cstdlib>
#include <vector>

struct Primitive {
    float f;
    Primitive() : f( drand48() ) {}
    bool intersect( Primitive *other )
    { return f * other->f > 0.6f; }
};

struct Box : public CGAL::Box_intersection_d::Unique_numbers
{
    Primitive *primitive;
    float lo[3], size;
    Box( Primitive *p ) : primitive( p ) {
        for( unsigned int d = 0; d < 3; ++d )
            lo[d] = 10.0 * drand48();
        size = 1 + 3.0 * drand48();
    }
};

struct Box_traits {
    typedef ::Box Box;
    typedef float Number_type;

    static Number_type get_lo( const Box& b, unsigned int dim )
    { return b.lo[ dim ]; }

    static Number_type get_hi( const Box& b, unsigned int dim )
    { return b.lo[ dim ] + b.size; }

    static unsigned int get_num( const Box& b )
    { return b.get_num();     }

    static unsigned int get_dim() { return 3; }
};

void fill_boxes( unsigned int n, std::vector<Box> &boxes ) {
    for( unsigned int i = 0; i < n; ++i )
        boxes.push_back( Box( new Primitive() ) );
}

void callback( const Box &a, const Box &b ) {
    if( a.primitive->intersect( b.primitive ) )
        std::cout << "intersection between box "
                  << a.get_num() << " and " << b.get_num() << std::endl;
};

int main() {
    std::vector<Box> boxes1, boxes2;
    fill_boxes( 100, boxes1 );
    fill_boxes( 100, boxes2 );
    box_intersection_d( boxes1.begin(), boxes1.end(),
                        boxes2.begin(), boxes2.end(), callback, Box_traits() );
}

