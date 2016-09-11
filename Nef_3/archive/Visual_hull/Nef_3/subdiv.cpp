// file: Visual_hull/Nef_3/subdiv.C

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <algorithm>
#include <vector>

using std::cout;
using std::cerr;
using std::cin;
using std::endl;

const double A = 1.0; // subdiv ratio along edge
const double B = 2.0;
const double Z = A+B;

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Plane_3                                      Plane;
typedef Kernel::Segment_3                                    Segment;
typedef Kernel::Tetrahedron_3                                Tetrahedron;

// vector of intermediate vertex coordinates
typedef std::vector< std::pair< Point, Point*> > Point_vector;


// custom items for Polyhedron, adds "cut" flag in vertex

template <class Refs, class T, class P>
class My_vertex : public CGAL::HalfedgeDS_vertex_base<Refs, T, P> {
public:
    bool cut;
    My_vertex() : cut( false) {} // repeat mandatory constructors
    My_vertex( const P& pt) 
        : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt), cut( false) {}
};

struct My_items : public CGAL::Polyhedron_items_3 {
    template <class Refs, class Traits>
    struct Vertex_wrapper {
        typedef typename Traits::Point_3  Point;
        typedef My_vertex<Refs, CGAL::Tag_true, Point> Vertex;
    };
};

typedef CGAL::Polyhedron_3<Kernel,My_items>                  Polyhedron;

typedef Polyhedron::Vertex                                   Vertex;
typedef Polyhedron::Vertex_iterator                          Vertex_iterator;
typedef Polyhedron::Halfedge_handle                          Halfedge_handle;
typedef Polyhedron::Edge_iterator                            Edge_iterator;
typedef Polyhedron::Facet_iterator                           Facet_iterator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;

void cut_degree_3( Polyhedron& P, Vertex_iterator vi) {
    Halfedge_handle h = vi->halfedge();
    // geometric precondition for a valid cut
    if ( CGAL::orientation( vi->point(), 
                            h->opposite()->vertex()->point(),
                            h->next()->vertex()->point(),
                            h->next()->opposite()->next()->vertex()->point()) 
        != CGAL::POSITIVE)
        return;
    // topological cut of the vertex
    P.split_vertex( h, h->next()->opposite());
    Halfedge_handle g = h->next()->opposite();
    P.split_vertex( g, g->next()->opposite());
    P.split_facet( h->next()->next()->opposite(), g->next());
    h->vertex()->cut = true;
    h->next()->vertex()->cut = true;
    g->next()->vertex()->cut = true;
}

double EPS = 0.000001;

void cutgeo_best_plane( Polyhedron& P, Point o, std::size_t n, 
                        Point_vector& new_points, std::size_t offset) {
    // Find geometrically optimal cut plane
    // Exhaustive search over all point triples to form plane
    std::size_t ci = 0, cj = 0, ck = 0;
    double sq_dist = -EPS;
    for ( std::size_t i = 0; i != n-2; ++i) {
        Point& p = new_points[ offset+i].first;
        for ( std::size_t j = i+1; j != n-1; ++j) {
            Point& q = new_points[ offset+j].first;
            for ( std::size_t k = j+1; k != n; ++k) {
                Point& r = new_points[ offset+k].first;
                std::size_t l = 0;
                for ( ; l != n; ++l) {
                    if ( l == i || l == j || l == k)
                        continue;
                    Tetrahedron t( new_points[ offset+l].first, p, q, r);
                    if ( t.volume() > EPS) // plane not separating o
                        break;
                }
                if ( l == n) { // p,q,r is a plane separating o from all others
                    Tetrahedron t( o, p, q, r);
                    if ( t.volume() > -EPS) { // should be true 
                        Plane pl( p, q, r);
                        double d = CGAL::squared_distance( pl, o);
                        if ( d > sq_dist) { // plane is better than previous
                            ci = i;
                            cj = j;
                            ck = k;
                            sq_dist = d;
                        }
                    }
                }
            }
        }
    }
    if ( cj == 0) { // Ooops, no plane found
        cerr << "Warning: couldn't resolve high-degree convex corner for "
            "numerical reasons." << endl;
        // hack: keep topological cut and don't smooth vertex coordinates
        // cleanup of new vertex vector
        new_points.erase( new_points.begin() + offset, new_points.end());
    } else {
        // prepare new coordinates
        // interpolate the three selected vertices
        Point& p = new_points[ offset+ci].first;
        Point& q = new_points[ offset+cj].first;
        Point& r = new_points[ offset+ck].first;
        Vector v = (o - CGAL::ORIGIN) * B;
        p = CGAL::ORIGIN + ((p - CGAL::ORIGIN) * A + v) / Z;
        q = CGAL::ORIGIN + ((q - CGAL::ORIGIN) * A + v) / Z;
        r = CGAL::ORIGIN + ((r - CGAL::ORIGIN) * A + v) / Z;
        // define plane and intersect all other edges with that plane
        // in case of expected numerical errors we set the vertex to o
        Plane pl( p, q, r);
        for ( std::size_t l = 0; l != n; ++l) {
            if ( l == ci || l == cj || l == ck)
                continue;
            Point& s = new_points[ offset+l].first;
            Segment seg( s, o);
            CGAL::Object rs = CGAL::intersection( seg, pl);
            if ( ! CGAL::assign( s, rs)) {
                cerr << "Warning: couldn't cut one edge of high-degree convex "
                    "corner properly." << endl;
                s = o; // fallback 
            }
        }
    }
}

void cutgeo_best_normal( Polyhedron& P, Point o, std::size_t n, 
                         Point_vector& new_points, std::size_t offset) {
    // Find geometrically optimal normal vector for cut plane
    // add negative inverse proportional locus vectors relative to o
    Vector v (0,0,0);
    for ( std::size_t i = 0; i != n; ++i) {
        Vector w = new_points[ offset+i].first - o;
        v = v + (w / (w*w));
    }
    // find point of smallest (positive) distance to o along vector v
    // if distance becomes negative, don't smooth new vertices
    std::size_t ci = 0;
    double dist = v * (new_points[ offset].first - o);
    for ( std::size_t i = 1; i != n; ++i) {
        double d = v * (new_points[ offset+i].first - o);
        if ( d < dist) {
            dist = d;
            ci = i;
        }
    }
    if ( dist < EPS) {
        cerr << "Warning: couldn't resolve high-degree convex corner for "
            "numerical reasons." << endl;
        // hack: keep topological cut and don't smooth vertex coordinates
        // cleanup of new vertex vector
        new_points.erase( new_points.begin() + offset, new_points.end());
        return;
    }
    // scale this shortest point
    Point q = CGAL::ORIGIN + ((new_points[ offset+ci].first - CGAL::ORIGIN) * A
               + ( o - CGAL::ORIGIN) * B) / Z;
    // define plane through q with normal v
    Plane pl( q, v);
    // assign new points as intersection with plane
    new_points[ offset+ci].first = q;
    for ( std::size_t l = 0; l != n; ++l) {
        if ( l == ci)
            continue;
        Point& s = new_points[ offset+l].first;
        Segment seg( s, o);
        CGAL::Object rs = CGAL::intersection( seg, pl);
        if ( ! CGAL::assign( s, rs)) {
            cerr << "Warning: couldn't cut one edge of high-degree convex "
                "corner properly." << endl;
            s = o; // fallback 
        }
    }
}

void cut_general( Polyhedron& P, Vertex_iterator vi, std::size_t n, 
                  Point_vector& new_points) {
    Point o = vi->point();
    Halfedge_handle h = vi->halfedge();
    // geometric precondition for a valid cut
    for ( std::size_t i = 0; i != n; ++i) {
        if ( CGAL::orientation( 
                 o,
                 h->opposite()->vertex()->point(),
                 h->next()->vertex()->point(),
                 h->next()->opposite()->next()->vertex()->point()) 
             != CGAL::POSITIVE) {
            Tetrahedron t( o,
                           h->opposite()->vertex()->point(),
                           h->next()->vertex()->point(),
                           h->next()->opposite()->next()->vertex()->point());
            CGAL_assertion( t.volume() <= 0);
            return;
        }
        h = h->next()->opposite();
    }
    std::size_t offset = new_points.size();
    // topological cut of the vertex
    Halfedge_handle hstart = h;
    h = h->next()->opposite();
    Halfedge_handle hfirst = h;
    while ( h != hstart) {
        Halfedge_handle g = h->next()->opposite();
        P.split_vertex( hstart, h);
        new_points.push_back( std::make_pair( h->opposite()->vertex()->point(),
                                              & (h->vertex()->point())));
        h = g;
    }
    new_points.push_back( std::make_pair( h->opposite()->vertex()->point(),
                                          & (h->vertex()->point())));
    CGAL_assertion( new_points.size() == offset + n);
    P.split_facet( h, hfirst->next()->opposite());
    //cutgeo_best_plane( P, o, n, new_points, offset);
    cutgeo_best_normal( P, o, n, new_points, offset);
}

void subdiv( Polyhedron& P) {
    if ( P.size_of_vertices() == 0)
        return;
    // Conservative size estimate assuming all vertices are convex.
    P.reserve( P.size_of_halfedges(), 
               2 * P.size_of_halfedges(), 
               P.size_of_facets() + P.size_of_vertices());
    // We use that new vertices/halfedges are appended at the end.
    Vertex_iterator last_v = P.vertices_end();
    -- last_v;  // the last of the old vertices
    Edge_iterator last_e = P.edges_end();
    -- last_e;  // the last of the old edges

    // store new coordinates for non-regular vertices in this vector
    Point_vector new_points;
    new_points.reserve(500);
    // first pass: cut regular vertices topologically
    Vertex_iterator vi = P.vertices_begin();
    do { 
        if ( vi->is_trivalent()) {
            cut_degree_3( P, vi);
        } else {
            std::size_t n = vi->vertex_degree();
            if ( n > 3)
                cut_general( P, vi, n, new_points);
        }
        if ( vi == last_v)
            break;
        ++vi;
    } while ( true);

    // second pass: assign new coordinates to "cut" regular vertices
    ++ last_e; // make it the past-the-end position again
    for ( Edge_iterator e = P.edges_begin(); e != last_e; ++e) {
        Vector v = e->vertex()->point() - CGAL::ORIGIN;
        Vector w = e->opposite()->vertex()->point() - CGAL::ORIGIN;
        if ( e->vertex()->cut)
            e->vertex()->point() = CGAL::ORIGIN + ((v*B+w*A)/Z);
        if ( e->opposite()->vertex()->cut)
            e->opposite()->vertex()->point() = CGAL::ORIGIN + ((v*A+w*B)/Z);
    };

    // third pass: assign new coordinates for non-regular vertices
    for ( Point_vector::iterator i = new_points.begin(); 
          i != new_points.end(); ++i) {
        * (i->second) = i->first;
    }
    CGAL_postcondition( P.is_valid());
}

void subdiv_regular( Polyhedron& P) {
    if ( P.size_of_vertices() == 0)
        return;
    // Conservative size estimate assuming all vertices are convex.
    P.reserve( P.size_of_halfedges(), 
               2 * P.size_of_halfedges(), 
               P.size_of_facets() + P.size_of_vertices());
    // We use that new vertices/halfedges are appended at the end.
    Vertex_iterator last_v = P.vertices_end();
    -- last_v;  // the last of the old vertices
    Edge_iterator last_e = P.edges_end();
    -- last_e;  // the last of the old edges

    // first pass: cut vertices topologically
    Vertex_iterator vi = P.vertices_begin();
    do { 
        if ( vi->is_trivalent())
            cut_degree_3( P, vi); // higher degree isn't convex here
        
        if ( vi == last_v)
            break;
        ++vi;
    } while ( true);
    
    // second pass: assign new coordinates to "cut" vertices
    ++ last_e; // make it the past-the-end position again
    for ( Edge_iterator e = P.edges_begin(); e != last_e; ++e) {
        Vector v = e->vertex()->point() - CGAL::ORIGIN;
        Vector w = e->opposite()->vertex()->point() - CGAL::ORIGIN;
        if ( e->vertex()->cut)
            e->vertex()->point() = CGAL::ORIGIN + ((v*B+w*A)/Z);
        if ( e->opposite()->vertex()->cut)
            e->opposite()->vertex()->point() = CGAL::ORIGIN + ((v*A+w*B)/Z);
    };
    CGAL_postcondition( P.is_valid());
}

int main( int argc, char* argv[]) {
    int iterations = 1;
    if ( argc > 1) {
        int i = atoi( argv[1]);
        if ( i > 1 && i < 100)
            iterations = i;
    }
    Polyhedron P;
    cin >> P;
    //P.normalize_border();
    //if ( P.size_of_border_edges() != 0) {
    //  cerr << "The input object has border edges. Cannot subdivide." << endl;
    //  std::exit(1);
    //}
    for ( int k = 1; k <= iterations; ++k) {
        cerr << k << ". iteration: ";
        CGAL::Timer t;
        t.start();
        if ( k == 1)
            subdiv( P);
        else
            subdiv_regular( P);
        t.stop();
        cerr << P.size_of_facets() << " facets in " 
             << t.time() << " seconds." << endl;
    }
    cout << P;
    return 0;
}
