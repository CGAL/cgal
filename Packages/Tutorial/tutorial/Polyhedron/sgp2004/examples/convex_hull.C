// file: examples/convex_hull.C

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Width_default_traits_3.h>
#include <CGAL/Width_3.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <algorithm>
#include <vector>


using std::cerr;
using std::endl;
using std::cout;
using std::cin;
using std::exit;

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;
typedef Polyhedron::Point_const_iterator                     Point_iterator;

typedef CGAL::Exact_predicates_exact_constructions_kernel    EKernel;
typedef CGAL::Polyhedron_3<EKernel>                          EPolyhedron;
typedef EKernel::Point_3                                     EPoint;
typedef CGAL::Width_default_traits_3<EKernel>                Width_traits;
typedef CGAL::Width_3<Width_traits>                          Width;


void convex_hull( const Polyhedron& P, Polyhedron& Q) {
    CGAL::convex_hull_3( P.points_begin(), P.points_end(), Q);
    cerr << "#vertices  : " << Q.size_of_vertices() << endl;
    cerr << "#facets    : " << Q.size_of_facets() << endl;
    cerr << "#edges     : " << (Q.size_of_halfedges() / 2) << endl;
}

void width( const Polyhedron& P) {
    std::vector< EPoint> epoints;
    for ( Point_iterator i = P.points_begin(); i != P.points_end(); ++i)
        epoints.push_back( EPoint( CGAL::to_double( i->x()),
                                   CGAL::to_double( i->y()),
                                   CGAL::to_double( i->z())));
    Width width( epoints.begin(), epoints.end());
    Width::RT num, denum;
    width.get_squared_width( num,denum);
    cerr << "width      : " << ( sqrt( CGAL::to_double( num) / 
                                       CGAL::to_double( denum))) << endl;
    cerr << "direction  : " << width.get_build_direction() << endl;
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

    user_time.reset();
    cerr << "Width ... " << endl;
    width( Q);
    cerr << "Width              : " << user_time.time() << " seconds." << endl;

    return 0;
}
