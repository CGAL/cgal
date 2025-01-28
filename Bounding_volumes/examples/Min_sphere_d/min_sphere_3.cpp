#include <CGAL/Simple_cartesian.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Random.h>

#include <iostream>
#include <cstdlib>

typedef CGAL::Simple_cartesian<double>           K;
typedef CGAL::Min_sphere_of_points_d_traits_3<K,double> Traits;
typedef CGAL::Min_sphere_of_spheres_d<Traits>    Min_sphere;
typedef K::Point_3                               Point;

const int n = 10;                        // number of points
const int d = 3;                         // dimension of points

int main ()
{
    Point         P[n];                  // n points
    CGAL::Random  r;                     // random number generator

    for (int i=0; i<n; ++i) {
        for (int j = 0; j < d; ++j) {
            P[i] = Point(r.get_double(), r.get_double(), r.get_double()); // random point
        }
    }

    Min_sphere  ms(P, P+n);             // smallest enclosing sphere

    Min_sphere::Cartesian_const_iterator ccib = ms.center_cartesian_begin(), ccie = ms.center_cartesian_end();
    std::cout << "center:";
    for( ; ccib != ccie; ++ccib){
      std::cout << " " << *ccib;
    }
    std::cout << std::endl << "radius: " << ms.radius() << std::endl;

    return 0;
}
