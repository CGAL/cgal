// examples/Polyhedron/polyhedron_prog_tetra_old.C
// -----------------------------------------------
#define CGAL_USE_POLYHEDRON_DESIGN_ONE 1
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>

typedef CGAL::Cartesian<double>      Kernel;
typedef Kernel::Point_3              Point;
typedef CGAL::Polyhedron_3<Kernel>   Polyhedron;
typedef Polyhedron::Vertex_iterator  Vertex_iterator;

int main() {
    Point p( 1.0, 0.0, 0.0);
    Point q( 0.0, 1.0, 0.0);
    Point r( 0.0, 0.0, 1.0);
    Point s( 0.0, 0.0, 0.0);

    Polyhedron P;
    P.make_tetrahedron( p, q, r, s);
    CGAL::set_ascii_mode( std::cout);
    for ( Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v)
        std::cout << v->point() << std::endl;
    return 0;
}
