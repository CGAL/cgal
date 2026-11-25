#include <CGAL/Simple_cartesian.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_2.h>
#include <CGAL/Random.h>

#include <iostream>
#include <array>

typedef  CGAL::Simple_cartesian<double>                   K;
typedef  CGAL::Min_sphere_of_points_d_traits_2<K,double>  Traits;
typedef  CGAL::Min_sphere_of_spheres_d<Traits>            Min_circle;
typedef  K::Point_2                                       Point;

int
main( int, char**)
{
    const int n = 100;
    std::array<Point, n> P;
    CGAL::Random  r;                     // random number generator

    for ( int i = 0; i < n; ++i){
      const double x = r.get_double();
      const double y = r.get_double();
      P.at(i) = Point(x, y);
    }

    Min_circle  mc( P.begin(), P.end());

    Min_circle::Cartesian_const_iterator ccib = mc.center_cartesian_begin(), ccie = mc.center_cartesian_end();
    std::cout << "center:";
    for( ; ccib != ccie; ++ccib){
      std::cout << " " << *ccib;
    }
    std::cout << std::endl << "radius: " << mc.radius() << std::endl;
    return 0;
}
