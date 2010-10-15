#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef Kernel::Point_3                    Point_3;
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;
typedef Polyhedron::Vertex_iterator        Vertex_iterator;

int main() {
    Point_3 p( 1.0, 0.0, 0.0);
    Point_3 q( 0.0, 1.0, 0.0);
    Point_3 r( 0.0, 0.0, 1.0);
    Point_3 s( 0.0, 0.0, 0.0);

    Polyhedron P;
    P.make_tetrahedron( p, q, r, s);
    CGAL::set_ascii_mode( std::cout);
    for ( Vertex_iterator v = P.vertices_begin(); v != P.vertices_end(); ++v)
        std::cout << v->point() << std::endl;
    return 0;
}
