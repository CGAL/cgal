#include <CGAL/Box_intersection_d.h>

#include <iostream>
#include <vector>
#include <CGAL/Random.h>

struct Primitive {
    float f;
    Primitive() : f( CGAL::default_random.get_double() ) {}
    bool intersect( Primitive *other )
    { return f * other->f > 0.6f; }
};

struct User_box : public CGAL::Box_intersection_d::Unique_numbers
{
    Primitive *primitive;
    float min[3], size;
    User_box( Primitive *p ) : primitive( p ) {
        for( unsigned int d = 0; d < 3; ++d )
            min[d] = 10.0 * CGAL::default_random.get_double();
        size = 1 + 3.0 * CGAL::default_random.get_double();
    }
};

struct Box_traits {
    typedef const User_box& Box;
    typedef float NT;

    static NT min( Box b, unsigned int dim )
    { return b.min[ dim ]; }

    static NT max( Box b, unsigned int dim )
    { return b.min[ dim ] + b.size; }

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

