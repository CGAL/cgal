// file: examples/mini_ball.C

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
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

//typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Cartesian< CGAL::Gmpq>                         Kernel;

//typedef CGAL::Simple_cartesian<double>                     Kernel;
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

void mini_ball( const Polyhedron& P) {
    Min_sphere min_sphere( P.points_begin(), P.points_end());
    cout << "Center point  : " << min_sphere.center() << endl;
    cout << "Square radius : " << min_sphere.squared_radius() << endl;
}

int main() {
    CGAL::Timer user_time;
    cerr << "Loading OFF file ... " << endl;
    user_time.start();
    Polyhedron P;
    cin >> P;
    cerr << "Loading OFF file   : " << user_time.time() << " seconds." << endl;

    user_time.reset();
    cerr << "Mini-ball ... " << endl;
    mini_ball( P);
    cerr << "Mini-ball          : " << user_time.time() << " seconds." << endl;
    return 0;
}
