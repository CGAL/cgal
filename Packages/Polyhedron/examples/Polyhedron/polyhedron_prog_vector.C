// polyhedron_prog_vector.C
// ------------------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/HalfedgeDS_using_vector.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>

typedef CGAL::Cartesian<double>                              R;
typedef CGAL::Polyhedron_default_traits_3<R>                 Traits;
typedef CGAL::Polyhedron_3< Traits, 
                            CGAL::Polyhedron_items_3, 
                            CGAL::HalfedgeDS_using_vector>   Polyhedron;
typedef Polyhedron::Point                                    Point;

int main() {
    Point p( 1.0, 0.0, 0.0);
    Point q( 0.0, 1.0, 0.0);
    Point r( 0.0, 0.0, 1.0);
    Point s( 0.0, 0.0, 0.0);

    Polyhedron P;    // alternative constructor: Polyhedron P(4,12,4);
    P.make_tetrahedron( p, q, r, s);
    CGAL::set_ascii_mode( std::cout);
    std::copy( P.points_begin(), P.points_end(),
	       std::ostream_iterator<Point>( std::cout, "\n"));
    return 0;
}
