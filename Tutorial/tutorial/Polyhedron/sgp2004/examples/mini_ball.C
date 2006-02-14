// file: examples/mini_ball.C

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_min_items_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <algorithm>
#include <vector>

#include <CGAL/Min_sphere_d.h>
#include <CGAL/Optimisation_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>



using std::cerr;
using std::endl;
using std::cout;
using std::cin;
using std::exit;

//typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
//typedef CGAL::Cartesian< CGAL::Gmpq>                         Kernel;

typedef CGAL::Simple_cartesian<double>                     Kernel;
typedef Kernel::Vector_3                                     Vector;
typedef Kernel::Point_3                                      Point;
//typedef CGAL::Polyhedron_3<Kernel>                           Polyhedron;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_min_items_3,
                           CGAL::HalfedgeDS_vector>        Polyhedron;

typedef Polyhedron::Point_const_iterator                     Point_iterator;

typedef CGAL::Optimisation_d_traits_3<Kernel>                 Traits;
typedef CGAL::Min_sphere_d<Traits>                            Min_sphere;
typedef CGAL::Min_sphere_of_spheres_d_traits_3<Kernel,double> TraitsS;
typedef CGAL::Min_sphere_of_spheres_d<TraitsS> Min_sphere_of_spheres;
typedef Min_sphere_of_spheres::Sphere                         Sphere;

void mini_ball( const Polyhedron& P) {
    Min_sphere min_sphere( P.points_begin(), P.points_end());
    Point p = min_sphere.center();
    cout << "Center point  : " 
         << CGAL::to_double(p.x()) << ' '
         << CGAL::to_double(p.y()) << ' '
         << CGAL::to_double(p.z()) << ' '
         << endl;
    cout << "Square radius : " << min_sphere.squared_radius() << endl;
    cout << "Double radius : " << sqrt( CGAL::to_double( 
                                    min_sphere.squared_radius())) << endl;
}

void min_sphere_of_spheres( const Polyhedron& P) {
    //std::vector<Sphere> S;
    Min_sphere_of_spheres min_sphere;
    min_sphere.prepare( P.size_of_vertices());
    //for ( Point_iterator i = P.points_begin(); i!= P.points_end(); ++i) {
    //    S.push_back( Sphere(*i, 0.0));
    //}
    for ( Point_iterator i = P.points_begin(); i!= P.points_end(); ++i) {
        min_sphere.insert( Sphere(*i, 0.0));
    }
    //Min_sphere_of_spheres min_sphere(S.begin(), S.end());
    Min_sphere_of_spheres::Support_iterator j = min_sphere.support_begin();
    int k = 0;
    while ( j != min_sphere.support_end()) {
        ++j;
        ++k;
    }
    cout << "#Spheres      : " << k << endl;
    cout << "Center point  : ";
    std::copy( min_sphere.center_cartesian_begin(), 
               min_sphere.center_cartesian_end(),
               std::ostream_iterator<double>( cout, " "));
    cout << endl;
    cout << "Double radius : " << min_sphere.radius() << endl;
}

int main() {
    CGAL::Timer user_time;
    cerr << "Loading OFF file ... " << endl;
    user_time.start();
    Polyhedron P;
    cin >> P;
    cerr << "Loading OFF file      : " << user_time.time()
         << " seconds." << endl;

    user_time.reset();
    cerr << "Mini-ball ... " << endl;
    mini_ball( P);
    cerr << "Mini-ball             : " << user_time.time() 
         << " seconds." << endl;

    user_time.reset();
    cerr << "Min_sphere_of_spheres ... " << endl;
    min_sphere_of_spheres( P);
    cerr << "Min_sphere_of_spheres : " << user_time.time()
         << " seconds." << endl;
    return 0;
}
