// polyhedron_prog_normals.C
// ------------------------------------------
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <iostream>
#include <algorithm>

struct Plane_equation {
    template <class Facet>
    void operator()( Facet& f) {
        typedef typename Facet::Plane            Plane;
        typedef typename Facet::Halfedge_handle  Halfedge_handle;
        Halfedge_handle h = f.halfedge();
        f.plane() = Plane( h->vertex()->point(),
                           h->next()->vertex()->point(),
                           h->next()->next()->vertex()->point());
    }
};

typedef CGAL::Homogeneous<int>                R;
typedef CGAL::Plane_3<R>                      Plane;
typedef CGAL::Polyhedron_default_traits_3<R>  Traits;
typedef CGAL::Polyhedron_3<Traits>            Polyhedron;
typedef Polyhedron::Point                     Point;

// The declaration of a plane iterator.
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

typedef Polyhedron::Facet                                     Facet;
typedef Polyhedron::Facet_iterator                            Facet_iterator;
typedef CGAL::Project_plane<Facet>                            Project_plane;
typedef CGAL::Iterator_project<Facet_iterator, Project_plane> Plane_iterator;

int main() {
    Point p( 1, 0, 0);
    Point q( 0, 1, 0);
    Point r( 0, 0, 1);
    Point s( 0, 0, 0);

    Polyhedron P;
    P.make_tetrahedron( p, q, r, s);
    std::for_each( P.facets_begin(), P.facets_end(), Plane_equation());

    CGAL::set_pretty_mode( std::cout);
    std::copy( Plane_iterator( P.facets_begin()), 
               Plane_iterator( P.facets_end()),
               std::ostream_iterator<Plane>( std::cout, "\n"));
    return 0;
}
