// file: examples/Min_ellipse_2/example_Min_ellipse_2.C

// includes
#include <CGAL/Homogeneous.h>
#include <CGAL/Point_2.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/Gmpz.h>

// typedefs
typedef  CGAL::Gmpz                       NT;
typedef  CGAL::Homogeneous<NT>            K;
typedef  CGAL::Point_2<K>                 Point;
typedef  CGAL::Min_ellipse_2_traits_2<K>  Traits;
typedef  CGAL::Min_ellipse_2<Traits>      Min_ellipse;

// main
int
main( int, char**)
{
    int     n = 100;
    Point*  P = new Point[ n];

    for ( int i = 0; i < n; ++i)
        P[ i] = Point( (i%2 == 0 ? i : -i), 0);
    // (0,0), (-1,0), (2,0), (-3,0), ...

    Min_ellipse  me1( P, P+n, false);    // very slow
    Min_ellipse  me2( P, P+n, true);     // fast

    CGAL::set_pretty_mode( std::cout);
    std::cout << me2;

    delete[] P;

    return( 0);
}

// ===== EOF ==================================================================
