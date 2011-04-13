// examples/Polyhedron/polyhedron_prog_normals_old.C
// -------------------------------------------------
#define CGAL_USE_POLYHEDRON_DESIGN_ONE 1
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <algorithm>

typedef CGAL::Homogeneous<int>       Kernel;
typedef Kernel::Point_3              Point;
typedef Kernel::Plane_3              Plane;
typedef CGAL::Polyhedron_3<Kernel>   Polyhedron;
typedef Polyhedron::Halfedge_handle  Halfedge_handle;
typedef Polyhedron::Facet            Facet;
typedef Polyhedron::Facet_iterator   Facet_iterator;

void compute_plane_equations( Facet& f) {
    Halfedge_handle h = f.halfedge();
    f.plane() = Plane( h->opposite()->vertex()->point(), 
		       h->vertex()->point(),
		       h->next()->vertex()->point());
};

int main() {
    Point p( 1, 0, 0);
    Point q( 0, 1, 0);
    Point r( 0, 0, 1);
    Point s( 0, 0, 0);

    Polyhedron P;
    P.make_tetrahedron( p, q, r, s);
    std::for_each( P.facets_begin(), P.facets_end(), compute_plane_equations);

    CGAL::set_pretty_mode( std::cout);
    Facet_iterator begin = P.facets_begin();
    for ( ; begin != P.facets_end(); ++begin)
        std::cout << begin->plane() << std::endl;
    return 0;
}
