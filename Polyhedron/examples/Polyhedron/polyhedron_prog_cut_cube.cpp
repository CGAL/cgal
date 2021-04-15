#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/polyhedron_cut_plane_3.h>
#include <iostream>

template <class Poly>
typename Poly::Halfedge_handle make_cube_3( Poly& P) {
    // appends a cube of size [0,1]^3 to the polyhedron P.
    CGAL_precondition( P.is_valid());
    typedef typename Poly::Point_3         Point;
    typedef typename Poly::Plane_3         Plane;
    typedef typename Poly::Halfedge_handle Halfedge_handle;
    Halfedge_handle h = P.make_tetrahedron( Point( 1, 0, 0),
                                            Point( 0, 0, 1),
                                            Point( 0, 0, 0),
                                            Point( 0, 1, 0));
    Halfedge_handle g = h->next()->opposite()->next();
    P.split_edge( h->next());
    P.split_edge( g->next());
    P.split_edge( g);
    h->next()->vertex()->point()     = Point( 1, 0, 1);
    g->next()->vertex()->point()     = Point( 0, 1, 1);
    g->opposite()->vertex()->point() = Point( 1, 1, 0);
    Halfedge_handle f = P.split_facet( g->next(), g->next()->next()->next());
    Halfedge_handle e = P.split_edge( f);
    e->vertex()->point() = Point( 1, 1, 1);
    P.split_facet( e, f->next()->next());
    CGAL_postcondition( P.is_valid());
    g = h;
    g->facet()->plane() = Plane( g->vertex()->point(),
                                 g->next()->vertex()->point(),
                                 g->next()->next()->vertex()->point());
    g = h->opposite();
    g->facet()->plane() = Plane( g->vertex()->point(),
                                 g->next()->vertex()->point(),
                                 g->next()->next()->vertex()->point());
    g = h->next()->opposite();
    g->facet()->plane() = Plane( g->vertex()->point(),
                                 g->next()->vertex()->point(),
                                 g->next()->next()->vertex()->point());
    g = h->next()->next()->opposite();
    g->facet()->plane() = Plane( g->vertex()->point(),
                                 g->next()->vertex()->point(),
                                 g->next()->next()->vertex()->point());
    g = h->next()->next()->next()->opposite();
    g->facet()->plane() = Plane( g->vertex()->point(),
                                 g->next()->vertex()->point(),
                                 g->next()->next()->vertex()->point());
    g = g->next()->next()->opposite();
    g->facet()->plane() = Plane( g->vertex()->point(),
                                 g->next()->vertex()->point(),
                                 g->next()->next()->vertex()->point());
    return h;
}

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Plane_3                                      Plane;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;
typedef Polyhedron::Halfedge_handle                          Halfedge_handle;

int main() {
    Polyhedron P;
    Halfedge_handle h = make_cube_3( P);
    Plane pl = Plane( Point( 0.5, 0.0, 0.0),
                      Point( 0.0, 0.0, 1.5),
                      Point( 0.0, 0.5, 0.0));
    CGAL::polyhedron_cut_plane_3( P, h, pl);
    std::cout << std::setprecision(17)<< P;
    return 0;
}
