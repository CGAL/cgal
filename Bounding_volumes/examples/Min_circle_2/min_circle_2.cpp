#include <CGAL/Simple_cartesian.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_2.h.h>

#include <iostream>

// typedefs
typedef  CGAL::Simple_cartesian<double>            K;
typedef  CGAL::Min_sphere_of_points_d_traits_2<K>  Traits;
typedef  CGAL::Min_sphere_of_spheres_d<Traits>     Min_circle;

typedef  K::Point_2                                Point;

int
main( int, char**)
{
    const int n = 100;
    Point P[n];

    for ( int i = 0; i < n; ++i){
      P[ i] = Point( (i%2 == 0 ? i : -i), 0, 1);
      // (0,0), (-1,0), (2,0), (-3,0), ...
    }

    Min_circle  mc( P, P+n);

    CGAL::set_pretty_mode( std::cout);
    std::cout << mc;

    return 0;
}
