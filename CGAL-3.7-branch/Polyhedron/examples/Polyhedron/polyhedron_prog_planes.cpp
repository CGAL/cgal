#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <algorithm>

struct Plane_equation {
    template <class Facet>
    typename Facet::Plane_3 operator()( Facet& f) {
        typename Facet::Halfedge_handle h = f.halfedge();
        typedef typename Facet::Plane_3  Plane;
        return Plane( h->vertex()->point(),
                      h->next()->vertex()->point(),
                      h->next()->next()->vertex()->point());
    }
};

typedef CGAL::Homogeneous<int>      Kernel;
typedef Kernel::Point_3             Point_3;
typedef Kernel::Plane_3             Plane_3;
typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;

int main() {
    Point_3 p( 1, 0, 0);
    Point_3 q( 0, 1, 0);
    Point_3 r( 0, 0, 1);
    Point_3 s( 0, 0, 0);
    Polyhedron P;
    P.make_tetrahedron( p, q, r, s);
    std::transform( P.facets_begin(), P.facets_end(), P.planes_begin(),
                    Plane_equation());
    CGAL::set_pretty_mode( std::cout);
    std::copy( P.planes_begin(), P.planes_end(),
               std::ostream_iterator<Plane_3>( std::cout, "\n"));
    return 0;
}
