// polyhedron_prog_normals.C
// -----------------------------------------------------------
#include <CGAL/Homogeneous.h>
#include <iostream>
#include <algorithm>
#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Homogeneous<int>                                R;
typedef CGAL::Halfedge_data_structure_polyhedron_default_3<R> HDS;
typedef CGAL::Polyhedron_default_traits_3<R>                  Traits;
typedef CGAL::Polyhedron_3<Traits,HDS>                        Polyhedron;
typedef Polyhedron::Point                                     Point;
typedef Polyhedron::Plane                                     Plane;
typedef Polyhedron::Halfedge_handle                           Halfedge_handle;
typedef Polyhedron::Facet                                     Facet;
typedef Polyhedron::Facet_iterator                            Facet_iterator;

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
