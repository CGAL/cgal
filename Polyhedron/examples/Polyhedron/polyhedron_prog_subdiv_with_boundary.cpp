#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <cctype>
#include <cmath>

typedef CGAL::Cartesian<double>                              Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;

typedef Polyhedron::Vertex                                   Vertex;
typedef Polyhedron::Vertex_iterator                          Vertex_iterator;
typedef Polyhedron::Halfedge_handle                          Halfedge_handle;
typedef Polyhedron::Edge_iterator                            Edge_iterator;
typedef Polyhedron::Facet_iterator                           Facet_iterator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;

void create_center_vertex( Polyhedron& P, Facet_iterator f) {
    Vector vec( 0.0, 0.0, 0.0);
    std::size_t order = 0;
    HF_circulator h = f->facet_begin();
    do {
        vec = vec + ( h->vertex()->point() - CGAL::ORIGIN);
        ++ order;
    } while ( ++h != f->facet_begin());
    CGAL_assertion( order >= 3); // guaranteed by definition of polyhedron
    Point center =  CGAL::ORIGIN + (vec / static_cast<double>(order));
    Halfedge_handle new_center = P.create_center_vertex( f->halfedge());
    new_center->vertex()->point() = center;
}

struct Smooth_old_vertex {
    Point operator()( const Vertex& v) const {
        std::size_t degree = CGAL::circulator_size( v.vertex_begin());
        if ( degree & 1) // odd degree only at border vertices
            return v.point();
        degree = degree / 2;
        double alpha = ( 4.0 - 2.0 * std::cos( 2.0 * CGAL_PI / degree)) / 9.0;
        Vector vec = (v.point() - CGAL::ORIGIN) * ( 1.0 - alpha);
        HV_circulator h = v.vertex_begin();
        do {
            // Even degree and border edges indicated non-manifold
            // vertex with two holes.
            if ( h->is_border_edge()) {
                std::cerr << "Error: non-manifold vertex. Erroneous smoothing."
                          << std::endl;
                return v.point();
            }
            vec = vec + ( h->opposite()->vertex()->point() - CGAL::ORIGIN)
              * alpha / static_cast<double>(degree);
            ++ h;
            CGAL_assertion( h != v.vertex_begin()); // even degree guaranteed
            ++ h;
        } while ( h != v.vertex_begin());
        return (CGAL::ORIGIN + vec);
    }
};

void flip_edge( Polyhedron& P, Halfedge_handle e) {
    if ( e->is_border_edge())
        return;
    Halfedge_handle h = e->next();
    P.join_facet( e);
    P.split_facet( h, h->next()->next());
}

void subdiv( Polyhedron& P) {
    if ( P.size_of_facets() == 0)
        return;
    // We use that new vertices/halfedges/facets are appended at the end.
    std::size_t nv = P.size_of_vertices();
    Vertex_iterator last_v = P.vertices_end();
    -- last_v;  // the last of the old vertices
    Edge_iterator last_e = P.edges_end();
    -- last_e;  // the last of the old edges
    Facet_iterator last_f = P.facets_end();
    -- last_f;  // the last of the old facets

    Facet_iterator f = P.facets_begin();    // create new center vertices
    do {
        create_center_vertex( P, f);
    } while ( f++ != last_f);

    std::vector<Point> pts;                    // smooth the old vertices
    pts.reserve( nv);  // get intermediate space for the new points
    ++ last_v; // make it the past-the-end position again
    std::transform( P.vertices_begin(), last_v, std::back_inserter( pts),
                    Smooth_old_vertex());
    std::copy( pts.begin(), pts.end(), P.points_begin());

    Edge_iterator e = P.edges_begin();              // flip the old edges
    ++ last_e; // make it the past-the-end position again
    while ( e != last_e) {
        Halfedge_handle h = e;
        ++e; // careful, incr. before flip since flip destroys current edge
        flip_edge( P, h);
    };
    CGAL_postcondition( P.is_valid());
}

void trisect_border_halfedge( Polyhedron& P, Halfedge_handle e) {
    CGAL_precondition( e->is_border());
    // Create two new vertices on e.
    e = e->prev();
    P.split_vertex( e, e->next()->opposite());
    P.split_vertex( e, e->next()->opposite());
    e = e->next();
    // We use later for the smoothing step that e->next()->next()
    // is our original halfedge we started with, i.e., its vertex is
    // from the unrefined mesh.  Split the face twice.
    Halfedge_handle h = e->opposite()->next();
    P.split_facet( e->next()->next()->opposite(), h);
    P.split_facet( e->next()->opposite(), h);
}

template <class OutputIterator>
void smooth_border_vertices( Halfedge_handle e, OutputIterator out) {
    CGAL_precondition( e->is_border());
    // We know that the vertex at this edge is from the unrefined mesh.
    // Get the locus vectors of the unrefined vertices in the neighborhood.
    Vector v0 = e->prev()->prev()->opposite()->vertex()->point() -CGAL::ORIGIN;
    Vector v1 = e->vertex()->point() - CGAL::ORIGIN;
    Vector v2 = e->next()->next()->next()->vertex()->point() - CGAL::ORIGIN;
    *out++ = CGAL::ORIGIN + (10.0 * v0 + 16.0 * v1 +        v2) / 27.0;
    *out++ = CGAL::ORIGIN + ( 4.0 * v0 + 19.0 * v1 +  4.0 * v2) / 27.0;
    *out++ = CGAL::ORIGIN + (       v0 + 16.0 * v1 + 10.0 * v2) / 27.0;
}

void subdiv_border( Polyhedron& P) {
    if ( P.size_of_facets() == 0)
        return;
    // We use that new halfedges are appended at the end.
    Edge_iterator last_e = P.edges_end();
    -- last_e;  // the last of the old edges
    Edge_iterator e = P.edges_begin(); // create trisected border edges
    do {
        if ( e->opposite()->is_border())
            trisect_border_halfedge( P, e->opposite());
        else if ( e->is_border())
            trisect_border_halfedge( P, e);
    } while ( e++ != last_e);
    e = P.edges_begin();               // smooth points on border edges
    std::vector<Point> pts;  // store new smoothed points temporarily
    do {
        if ( e->opposite()->is_border())
            smooth_border_vertices( e->opposite(), std::back_inserter(pts));
        else if ( e->is_border())
            smooth_border_vertices( e, std::back_inserter(pts));
    } while ( e++ != last_e);
    e = P.edges_begin(); // copy smoothed points back
    std::vector<Point>::iterator i = pts.begin();
    do {
        if ( e->opposite()->is_border()) {
            e->vertex()->point() = *i++;
            e->opposite()->vertex()->point() = *i++;
            e->opposite()->next()->vertex()->point() = *i++;
        } else if ( e->is_border()) {
            e->opposite()->vertex()->point() = *i++;
            e->vertex()->point() = *i++;
            e->next()->vertex()->point() = *i++;
        }
    } while ( e++ != last_e);
    CGAL_assertion( i == pts.end());
    CGAL_postcondition( P.is_valid());
}

using namespace std;

int main( int argc, char* argv[]) {
    if ( argc > 2 || (argc == 2 && ! isdigit( argv[1][0]))) {
        cerr << "Usage: " << argv[0] << " [<n>]" << endl;
        cerr << "    subdivides <n> times the polyhedron read from stdin."
             << endl;
        exit(1);
    }
    int n = 1;
    if ( argc >= 2)
        n = atoi( argv[1]);
    if ( n < 1 || n > 12) {
        cerr << "Error: Choose reasonable value for <n> in [1..12]" << endl;
        exit(1);
    }
    Polyhedron P;
    cin >> P;
    for ( int i = 0; i != n; ++i) {
        cerr << "Subdivision " << i+1 << " ..." << endl;
        subdiv( P);
        if ( i & 1)
            subdiv_border( P);
    }
    cout << P;
    return 0;
}
