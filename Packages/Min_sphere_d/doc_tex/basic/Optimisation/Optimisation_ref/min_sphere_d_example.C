
#include<CGAL/Cartesian.h>
#include<iostream>
#include<cstdlib>
#include<CGAL/Random.h>
#include<CGAL/Optimisation_d_traits_d.h>
#include<CGAL/Min_sphere_d.h>

typedef CGAL::Cartesian<double>                R;
typedef CGAL::Optimisation_d_traits_d<R>       Traits;
typedef CGAL::Min_sphere_d<Traits>             Min_sphere;
typedef R::Point_d                             Point;

const int n = 10;                        // number of points
const int d = 5;                         // dimension of points

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

