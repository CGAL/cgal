// examples/Polyhedron/polyhedron_prog_point_iterator_old.C
// --------------------------------------------------------
#define CGAL_USE_POLYHEDRON_DESIGN_ONE 1
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <iterator>
#include <algorithm>

typedef CGAL::Cartesian<double>                              Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;

// The declaration of a point iterator for a given Polyhedron
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

typedef Polyhedron::Vertex                                Vertex;
typedef Polyhedron::Vertex_iterator                       Vertex_iterator;
typedef CGAL::Project_point<Vertex>                       Project;
typedef CGAL::Iterator_project<Vertex_iterator, Project>  Point_iterator;

int main() {
    Point p( 0.0, 0.0, 0.0);
    Point q( 1.0, 0.0, 0.0);
    Point r( 0.0, 1.0, 0.0);
    Point s( 0.0, 0.0, 1.0);

    Polyhedron P;
    P.make_tetrahedron( p, q, r, s);
    CGAL::set_ascii_mode( std::cout);
    // explicit use in a loop
    Point_iterator begin = P.vertices_begin();
    for ( ; begin != P.vertices_end(); ++begin)
        std::cout << "(" << (*begin) << ") ";
    std::cout << std::endl;
    // used with a generic algorithm
    std::copy( Point_iterator( P.vertices_begin()),
	       Point_iterator( P.vertices_end()),
	       std::ostream_iterator<Point>( std::cout, ",  "));
    std::cout << std::endl;
    return 0;
}
