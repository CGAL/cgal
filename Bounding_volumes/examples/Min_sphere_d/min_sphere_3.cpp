#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <cstdlib>
#include <CGAL/Random.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>

typedef CGAL::Simple_cartesian<double>           K;
typedef CGAL::Min_sphere_of_points_d_traits_3<K> Traits;
typedef CGAL::Min_sphere_of_spheres_d<Traits>    Min_sphere;
typedef K::Point_3                               Point;

const int n = 10;                        // number of points
const int d = 3;                         // dimension of points

int main ()
{
    Point         P[n];                  // n points
    double        coord[d];              // d coordinates
    CGAL::Random  r;                     // random number generator

    for (int i=0; i<n; ++i) {
        for (int j=0; j<d; ++j)
            coord[j] = r.get_double();
        P[i] = Point(d, coord, coord+d); // random point
    }

    Min_sphere  ms (P, P+n);             // smallest enclosing sphere

    CGAL::set_pretty_mode (std::cout);
    std::cout << ms;                     // output the sphere

    return 0;
}
