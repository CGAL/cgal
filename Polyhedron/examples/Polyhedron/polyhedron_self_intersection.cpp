#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>

using std::cerr;
using std::endl;
using std::cout;
using std::cin;
using std::exit;

typedef CGAL::Bbox_3                                     Bbox;
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Triangle_3                               Triangle;
typedef Kernel::Segment_3                                Segment;
typedef CGAL::Polyhedron_3<Kernel>                       Polyhedron;
typedef Polyhedron::Halfedge_const_handle                Halfedge_const_handle;
typedef Polyhedron::Facet_const_iterator                 Facet_const_iterator;
typedef Polyhedron::Facet_const_handle                   Facet_const_handle;
typedef CGAL::Box_intersection_d::Box_with_handle_d<
            double, 3, Facet_const_handle>               Box;

std::vector<Triangle> triangles;

struct Intersect_facets {
    void operator()( const Box* b, const Box* c) const {
        Halfedge_const_handle h = b->handle()->halfedge();
        // check for shared egde --> no intersection
        if ( h->opposite()->facet() == c->handle()
             || h->next()->opposite()->facet() == c->handle()
             || h->next()->next()->opposite()->facet() == c->handle())
            return;
        // check for shared vertex --> maybe intersection, maybe not
        Halfedge_const_handle g = c->handle()->halfedge();
        Halfedge_const_handle v;
        if ( h->vertex() == g->vertex())
            v = g;
        if ( h->vertex() == g->next()->vertex())
            v = g->next();
        if ( h->vertex() == g->next()->next()->vertex())
            v = g->next()->next();
        if ( v == Halfedge_const_handle()) {
            h = h->next();
            if ( h->vertex() == g->vertex())
                v = g;
            if ( h->vertex() == g->next()->vertex())
                v = g->next();
            if ( h->vertex() == g->next()->next()->vertex())
                v = g->next()->next();
            if ( v == Halfedge_const_handle()) {
                h = h->next();
                if ( h->vertex() == g->vertex())
                    v = g;
                if ( h->vertex() == g->next()->vertex())
                    v = g->next();
                if ( h->vertex() == g->next()->next()->vertex())
                    v = g->next()->next();
            }
        }
        if ( v != Halfedge_const_handle()) {
            // found shared vertex:
            CGAL_assertion( h->vertex() == v->vertex());
            // geomtric check if the opposite segments intersect the triangles
            Triangle t1( h->vertex()->point(),
                         h->next()->vertex()->point(),
                         h->next()->next()->vertex()->point());
            Triangle t2( v->vertex()->point(),
                         v->next()->vertex()->point(),
                         v->next()->next()->vertex()->point());
            Segment  s1( h->next()->vertex()->point(),
                         h->next()->next()->vertex()->point());
            Segment  s2( v->next()->vertex()->point(),
                         v->next()->next()->vertex()->point());
            if ( CGAL::do_intersect( t1, s2)) {
                //cerr << "Triangles intersect (t1,s2):\n    T1: " << t1
                //     << "\n    T2 :" << t2 << endl;
                triangles.push_back(t1);
                triangles.push_back(t2);
            } else if ( CGAL::do_intersect( t2, s1)) {
                //cerr << "Triangles intersect (t2,s1):\n    T1: " << t1
                //     << "\n    T2 :" << t2 << endl;
                triangles.push_back(t1);
                triangles.push_back(t2);
            }
            return;
        }
        // check for geometric intersection
        Triangle t1( h->vertex()->point(),
                     h->next()->vertex()->point(),
                     h->next()->next()->vertex()->point());
        Triangle t2( g->vertex()->point(),
                     g->next()->vertex()->point(),
                     g->next()->next()->vertex()->point());
        if ( CGAL::do_intersect( t1, t2)) {
            //cerr << "Triangles intersect:\n    T1: " << t1 << "\n    T2 :"
            //     << t2 << endl;
            triangles.push_back(t1);
            triangles.push_back(t2);
        }
    }
};

void write_off() {
    cout << "OFF\n" << (triangles.size() * 3) << ' ' << triangles.size()
         << " 0\n";
    for ( std::vector<Triangle>::iterator i = triangles.begin();
          i != triangles.end(); ++i) {
        cout << i->vertex(0) << '\n';
        cout << i->vertex(1) << '\n';
        cout << i->vertex(2) << '\n';
    }
    for ( std::size_t k = 0; k != triangles.size(); ++k) {
        cout << "3  " << (3*k) << ' ' << (3*k+1) << ' ' << (3*k+2) << '\n';
    }
}

void intersection( const Polyhedron& P) {
    std::vector<Box> boxes;
    boxes.reserve( P.size_of_facets());
    for ( Facet_const_iterator i = P.facets_begin(); i != P.facets_end(); ++i){
        boxes.push_back(
            Box( i->halfedge()->vertex()->point().bbox()
               + i->halfedge()->next()->vertex()->point().bbox()
               + i->halfedge()->next()->next()->vertex()->point().bbox(),
                 i));
    }
    std::vector<const Box*> box_ptr;
    box_ptr.reserve( P.size_of_facets());
    for ( std::vector<Box>::iterator j = boxes.begin(); j != boxes.end(); ++j){
        box_ptr.push_back( &*j);
    }
    CGAL::box_self_intersection_d( box_ptr.begin(), box_ptr.end(),
                                   Intersect_facets(), std::ptrdiff_t(2000));
}

int main(int argc, char* argv[]) {
    CGAL::Timer user_time;
    cerr << "Loading OFF file ... " << endl;
    user_time.start();
    Polyhedron P;
    std::ifstream in1((argc>1)?argv[1]:"data/tetra_intersected_by_triangle.off");
    in1 >> P;
    cerr << "Loading OFF file   : " << user_time.time() << " seconds." << endl;
    if ( ! P.is_pure_triangle()) {
        cerr << "The input object is not triangulated. Cannot intersect."
                  << endl;
        exit(1);
    }
    user_time.reset();
    cerr << "Intersection ... " << endl;
    intersection( P);
    cerr << "Intersection       : " << user_time.time() << " seconds." << endl;
    write_off();

    return 0;
}
