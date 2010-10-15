#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double>               Kernel;
typedef Kernel::Point_3                              Point_3;
typedef CGAL::Polyhedron_3<Kernel>                   Polyhedron;
typedef Polyhedron::Facet_iterator                   Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

int main() {
    Point_3 p( 0.0, 0.0, 0.0);
    Point_3 q( 1.0, 0.0, 0.0);
    Point_3 r( 0.0, 1.0, 0.0);
    Point_3 s( 0.0, 0.0, 1.0);

    Polyhedron P;
    P.make_tetrahedron( p, q, r, s);

    // Write polyhedron in Object File Format (OFF).
    CGAL::set_ascii_mode( std::cout);
    std::cout << "OFF" << std::endl << P.size_of_vertices() << ' '
              << P.size_of_facets() << " 0" << std::endl;
    std::copy( P.points_begin(), P.points_end(),
               std::ostream_iterator<Point_3>( std::cout, "\n"));
    for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
        Halfedge_facet_circulator j = i->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        std::cout << CGAL::circulator_size(j) << ' ';
        do {
            std::cout << ' ' << std::distance(P.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
        std::cout << std::endl;
    }
    return 0;
}
