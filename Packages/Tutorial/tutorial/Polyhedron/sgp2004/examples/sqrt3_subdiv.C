// file: examples/sqrt3_subdiv.C
// Revised and simplified to work with CGAL 3.1 and triangulated meshes only.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Real_timer.h>
#include <iostream>
#include <algorithm>
#include <vector>

using std::cerr;
using std::endl;
using std::cout;
using std::cin;
using std::exit;

// no plane equations in facets
class Polyhedron_min_items_3 {
public:
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef typename Traits::Point_3 Point;
        typedef CGAL::HalfedgeDS_vertex_base< Refs, CGAL::Tag_true, Point> 
                                                                   Vertex;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper {
        // with prev pointer, 1.05 s, without 0.94
        typedef CGAL::HalfedgeDS_halfedge_base< Refs, CGAL::Tag_false>
                                                                   Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        // with plane and prev pointer. 1.20 s
        typedef typename Traits::Plane_3 Plane;
        typedef CGAL::HalfedgeDS_face_base< Refs, CGAL::Tag_true>  Face;
    };
};


typedef CGAL::Simple_cartesian<float>                       Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel, Polyhedron_min_items_3,
                           CGAL::HalfedgeDS_vector>   Polyhedron;

typedef Polyhedron::Vertex                                   Vertex;
typedef Polyhedron::Vertex_iterator                          Vertex_iterator;
typedef Polyhedron::Halfedge_handle                          Halfedge_handle;
typedef Polyhedron::Edge_iterator                            Edge_iterator;
typedef Polyhedron::Facet_iterator                           Facet_iterator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;

void create_centroid( Polyhedron& P, Facet_iterator f) {
    Halfedge_handle h = f->halfedge();
    Vector vec = h->vertex()->point() - CGAL::ORIGIN;
    vec = vec + (h->next()->vertex()->point() - CGAL::ORIGIN);
    vec = vec + (h->next()->next()->vertex()->point() - CGAL::ORIGIN);
    Halfedge_handle new_center = P.create_center_vertex( h);
    new_center->vertex()->point() = CGAL::ORIGIN + (vec / 3.0);
}

struct Smooth_old_vertex {
    Point operator()( const Vertex& v) const {
        CGAL_precondition((CGAL::circulator_size( v.vertex_begin()) & 1) == 0);
        std::size_t degree = CGAL::circulator_size( v.vertex_begin()) / 2;
        double alpha = ( 4.0 - 2.0 * cos( 2.0 * CGAL_PI / degree)) / 9.0;
        Vector vec = (v.point() - CGAL::ORIGIN) * ( 1.0 - alpha);
        HV_circulator h = v.vertex_begin();
        do {
            vec = vec + ( h->opposite()->vertex()->point() - CGAL::ORIGIN) 
                       * alpha / degree;
            ++ h;
            CGAL_assertion( h != v.vertex_begin()); // even degree guaranteed
            ++ h;
        } while ( h != v.vertex_begin());
        return (CGAL::ORIGIN + vec);
    }
};

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
        create_centroid( P, f);
    } while ( f++ != last_f);

    std::vector<Point> pts;                    // smooth the old vertices
    pts.reserve( nv);  // get intermediate space for the new points
    ++ last_v; // make it the past-the-end position again
    std::transform( P.vertices_begin(), last_v, std::back_inserter( pts), 
                    Smooth_old_vertex());
    std::copy( pts.begin(), pts.end(), P.points_begin());

    ++ last_e; // make it the past-the-end position again
    for ( Edge_iterator e = P.edges_begin(); e != last_e; ++e)
        P.flip_edge(e);    // flip the old edges
    CGAL_postcondition( P.is_valid());
}

int main() {
    CGAL::Real_timer runtime;
    cerr << "Loading OFF file ... " << endl;
    runtime.start();
    Polyhedron P;
    cin >> P;
    cerr << "Loading OFF file   : " << runtime.time() << " seconds." << endl;

    P.normalize_border();
    if ( P.size_of_border_edges() != 0) {
        cerr << "The input object has border edges. Cannot subdivide." 
                  << endl;
        exit(1);
    }
    if ( ! P.is_pure_triangle()) {
        cerr << "The input object is not triangulated. Cannot subdivide." 
                  << endl;
        exit(1);
    }
    P.reserve( P.size_of_vertices() + P.size_of_facets(),
               P.size_of_halfedges() + 6 * P.size_of_facets(),
               3 * P.size_of_facets());
    runtime.reset();
    cerr << "Sqrt-3 Subdiv ... " << endl;
    subdiv( P);
    cerr << "Sqrt-3 Subdiv      : " << runtime.time() << " seconds." << endl;

    runtime.reset();
    cerr << "Saving OFF file ... " << endl;
    cout << P;
    cerr << "Saving OFF file    : " << runtime.time() << " seconds." << endl;
    return 0;
}
