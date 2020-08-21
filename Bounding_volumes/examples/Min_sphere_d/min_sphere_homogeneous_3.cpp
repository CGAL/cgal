#include <CGAL/Exact_integer.h>
#include <CGAL/Homogenous.h>
#include <CGAL/Random.h>
#include <CGAL/Min_sphere_annulus_d_traits_d.h>
#include <CGAL/Min_sphere_d.h>

#include <iostream>
#include <cstdlib>

typedef CGAL::Homogeneous<CGAL::Exact_integer>   K;
typedef CGAL::Min_sphere_annulus_d_traits_d<K>   Traits;
typedef CGAL::Min_sphere_d<Traits>               Min_sphere;
typedef K::Point_3                               Point;

const int n = 10;                        // number of points

int main ()
{
    Point         P[n];                  // n points
    CGAL::Random  r;                     // random number generator

    for (int i=0; i<n; ++i) {
      P[i] = Point(r.get_int(0, 1000),r.get_int(0, 1000), r.get_int(0, 1000), 1 );
    }

    Min_sphere  ms (P, P+n);             // smallest enclosing sphere

    CGAL::set_pretty_mode (std::cout);
    std::cout << ms;                     // output the sphere

    return 0;
}
