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

struct User_box : public CGAL::Box_intersection_d::Unique_numbers
{
    Primitive *primitive;
    float lo[3], size;
    User_box( Primitive *p ) : primitive( p ) {
        for( unsigned int d = 0; d < 3; ++d )
            lo[d] = 10.0 * drand48();
        size = 1 + 3.0 * drand48();
    }
};

struct Box_traits {
    typedef const User_box& Box;
    typedef float NT;

    static NT get_lo( Box b, unsigned int dim )
    { return b.lo[ dim ]; }

    static NT get_hi( Box b, unsigned int dim )
    { return b.lo[ dim ] + b.size; }

    static unsigned int get_id( Box b )
    { return b.get_id();     }

    static unsigned int get_dim() { return 3; }
};

void fill_boxes( unsigned int n, std::vector<User_box> &boxes ) {
    for( unsigned int i = 0; i < n; ++i )
        boxes.push_back( User_box( new Primitive() ) );
}

void callback( const User_box &a, const User_box &b ) {
    if( a.primitive->intersect( b.primitive ) )
        std::cout << "intersection between box "
                  << a.get_id() << " and " << b.get_id() << std::endl;
};

int main() {
    std::vector<User_box> a, b;
    fill_boxes( 100, a );
    fill_boxes( 100, b );
    CGAL::box_intersection_d_custom_traits(
          a.begin(), a.end(),
          b.begin(), b.end(), callback, Box_traits() );
}

