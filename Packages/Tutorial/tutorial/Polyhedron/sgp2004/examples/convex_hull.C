// file: examples/convex_hull.C

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <algorithm>
#include <vector>

#include <CGAL/Min_sphere_d.h>
#include <CGAL/Optimisation_d_traits_3.h>


using std::cerr;
using std::endl;
using std::cout;
using std::cin;
using std::exit;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
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

typedef CGAL::Optimisation_d_traits_3<Kernel>                Traits;
typedef CGAL::Min_sphere_d<Traits>                           Min_sphere;

void convex_hull( const Polyhedron& P, Polyhedron& Q) {
    CGAL::convex_hull_3( P.points_begin(), P.points_end(), Q);
    cerr << "#vertices  : " << Q.size_of_vertices() << endl;
    cerr << "#facets    : " << Q.size_of_facets() << endl;
    cerr << "#edges     : " << (Q.size_of_halfedges() / 2) << endl;
}

int main() {
    CGAL::Timer user_time;
    cerr << "Loading OFF file ... " << endl;
    user_time.start();
    Polyhedron P;
    cin >> P;
    cerr << "Loading OFF file   : " << user_time.time() << " seconds." << endl;

    user_time.reset();
    cerr << "Convex hull ... " << endl;
    Polyhedron Q;
    convex_hull( P, Q);
    cerr << "Convex hull        : " << user_time.time() << " seconds." << endl;

    user_time.reset();
    cerr << "Saving OFF file ... " << endl;
    cout << Q;
    cerr << "Saving OFF file    : " << user_time.time() << " seconds." << endl;
    return 0;
}
