// examples/Polyhedron/polyhedron_prog_off_old.C
// ---------------------------------------------
#define CGAL_USE_POLYHEDRON_DESIGN_ONE 1
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>
#include <iostream>

typedef CGAL::Cartesian<double>                          Kernel;
typedef Kernel::Point_3                                  Point;
typedef CGAL::Polyhedron_3<Kernel>                       Polyhedron;
typedef Polyhedron::Vertex                               Vertex;
typedef Polyhedron::Vertex_iterator                      Vertex_iterator;
typedef Polyhedron::Facet_iterator                       Facet_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator
                                            Halfedge_around_facet_circulator;

typedef CGAL::Project_point<Vertex>                           Project;
typedef CGAL::Iterator_project<Vertex_iterator, Project>      Point_iterator;

int main() {
    Point p( 0.0, 0.0, 0.0);
    Point q( 1.0, 0.0, 0.0);
    Point r( 0.0, 1.0, 0.0);
    Point s( 0.0, 0.0, 1.0);

    Polyhedron P;
    P.make_tetrahedron( p, q, r, s);

    // Write polyhedron on Object File Format (OFF).
    CGAL::set_ascii_mode( std::cout);
    std::cout << "OFF" << std::endl;
    std::cout << P.size_of_vertices() << ' ' 
	      << P.size_of_facets() << " 0" << std::endl;
    std::copy( Point_iterator( P.vertices_begin()), 
	       Point_iterator( P.vertices_end()), 
	       std::ostream_iterator<Point>( std::cout, "\n"));
    for ( Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
        Halfedge_around_facet_circulator j = i->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        std::cout << CGAL::circulator_size(j) << " ";
        do {
            std::cout << " " << std::distance(P.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
        std::cout << std::endl;
    }
    return 0;
}
