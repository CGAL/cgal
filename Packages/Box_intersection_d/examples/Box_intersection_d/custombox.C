#include <CGAL/box_intersection_d.h>

#include <iostream>
#include <vector>
#include <CGAL/Random.h>

struct Primitive {
    double f;
    Primitive() : f( CGAL::default_random.get_double() ) {}
    bool intersect( Primitive *other )
    { return f * other->f > 0.6; }
};

struct Box : public CGAL::Box_intersection_d::Box_d< double, 3 >
{
    typedef CGAL::Box_intersection_d::Box_d< double, 3 > Base;
    Primitive *primitive;
    Box( Primitive *p ) : primitive( p ) {
        double d1 = 10.0 * CGAL::default_random.get_double();
        double d2 = 10.0 * CGAL::default_random.get_double();
        double d3 = 10.0 * CGAL::default_random.get_double();
        *((Base*)(this)) = Base( 
            CGAL::Bbox_3( d1, d2, d3,
                          d1 + 1.0 + CGAL::default_random.get_double(),
                          d2 + 1.0 + CGAL::default_random.get_double(),
                          d3 + 1.0 + CGAL::default_random.get_double()));
    }
};

void fill_boxes( unsigned int n, std::vector< Box >& boxes ) {
    for( unsigned int i = 0; i < n; ++i )
        boxes.push_back( Box( new Primitive() ) );
}

void callback( const Box &a, const Box &b ) {
    if( a.primitive->intersect( b.primitive ) )
        std::cout << "intersection between box "
            << a.id() << " and " << b.id() << std::endl;
};

int main() {
    std::vector< Box > boxes1, boxes2;
    fill_boxes( 100, boxes1 );
    fill_boxes( 100, boxes2 );
    CGAL::box_intersection_d( boxes1.begin(), boxes1.end(),
                              boxes2.begin(), boxes2.end(), callback );
}

