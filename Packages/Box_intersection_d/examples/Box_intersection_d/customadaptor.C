#include <CGAL/box_intersection_d.h>

#include <iostream>
#include <vector>
#include <CGAL/Random.h>

struct Primitive {
    float f;
    Primitive() : f( CGAL::default_random.get_double() ) {}
    bool intersect( Primitive *other )
    { return f * other->f > 0.6f; }
};

struct User_box : public CGAL::Box_intersection_d::Unique_numbers<int>
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
    typedef float           NT;
    typedef std::size_t     Id_type;

    static NT min_coord( Box b, int dim )  { return b.min[ dim ]; }
    static NT max_coord( Box b, int dim )  { return b.min[ dim ] + b.size; }
    static std::size_t id( Box b )         { return b.id(); }
    static const int dimension()           { return 3; }
};

void fill_boxes( unsigned int n, std::vector<User_box> &boxes ) {
    for( unsigned int i = 0; i < n; ++i )
        boxes.push_back( User_box( new Primitive() ) );
}

void callback( const User_box &a, const User_box &b ) {
    if( a.primitive->intersect( b.primitive ) )
        std::cout << "intersection between box "
                  << a.id() << " and " << b.id() << std::endl;
};

int main() {
    std::vector<User_box> a, b;
    fill_boxes( 100, a );
    fill_boxes( 100, b );
    CGAL::box_intersection_d( a.begin(), a.end(), b.begin(), b.end(),
                              callback, Box_traits() );
}

