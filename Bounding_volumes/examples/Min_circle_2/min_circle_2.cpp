#include <CGAL/Simple_cartesian.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_2.h>
#include <CGAL/Random.h>

#include <iostream>

typedef  CGAL::Simple_cartesian<double>                   K;
typedef  CGAL::Min_sphere_of_points_d_traits_2<K,double>  Traits;
typedef  CGAL::Min_sphere_of_spheres_d<Traits>            Min_circle;
typedef  K::Point_2                                       Point;

int
main( int, char**)
{
    const int n = 100;
    Point P[n];
    CGAL::Random  r;                     // random number generator

    for ( int i = 0; i < n; ++i){
      P[ i] = Point(r.get_double(), r.get_double());
    }

    Min_circle  mc( P, P+n);

    Min_circle::Cartesian_const_iterator ccib = mc.center_cartesian_begin(), ccie = mc.center_cartesian_end();
    std::cout << "center:";
    for( ; ccib != ccie; ++ccib){
      std::cout << " " << *ccib;
    }
    std::cout << std::endl << "radius: " << mc.radius() << std::endl;
    return 0;
}
