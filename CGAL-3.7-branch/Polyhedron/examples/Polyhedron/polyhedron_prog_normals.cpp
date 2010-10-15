#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <algorithm>

struct Normal_vector {
    template <class Facet>
    typename Facet::Plane_3 operator()( Facet& f) {
        typename Facet::Halfedge_handle h = f.halfedge();
        // Facet::Plane_3 is the normal vector type. We assume the
        // CGAL Kernel here and use its global functions.
        return CGAL::cross_product(
          h->next()->vertex()->point() - h->vertex()->point(),
          h->next()->next()->vertex()->point() - h->next()->vertex()->point());
    }
};

typedef CGAL::Homogeneous<int>                         Kernel;
typedef Kernel::Point_3                                Point_3;
typedef Kernel::Vector_3                               Vector_3;
typedef CGAL::Polyhedron_traits_with_normals_3<Kernel> Traits;
typedef CGAL::Polyhedron_3<Traits>                     Polyhedron;

int main() {
    Point_3 p( 1, 0, 0);
    Point_3 q( 0, 1, 0);
    Point_3 r( 0, 0, 1);
    Point_3 s( 0, 0, 0);
    Polyhedron P;
    P.make_tetrahedron( p, q, r, s);
    std::transform( P.facets_begin(), P.facets_end(), P.planes_begin(),
                    Normal_vector());
    CGAL::set_pretty_mode( std::cout);
    std::copy( P.planes_begin(), P.planes_end(),
               std::ostream_iterator<Vector_3>( std::cout, "\n"));
    return 0;
}
