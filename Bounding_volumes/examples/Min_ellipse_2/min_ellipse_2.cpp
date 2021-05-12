#include <CGAL/Cartesian.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/Exact_rational.h>

#include <cassert>

typedef  CGAL::Exact_rational             NT;
typedef  CGAL::Cartesian<NT>              K;
typedef  CGAL::Point_2<K>                 Point;
typedef  CGAL::Min_ellipse_2_traits_2<K>  Traits;
typedef  CGAL::Min_ellipse_2<Traits>      Min_ellipse;


int
main( int, char**)
{
    const int n = 200;
    Point     P[n];

    for ( int i = 0; i < n; ++i)
        P[ i] = Point( i % 2 ? i : -i , 0);
    // (0,0), (-1,0), (2,0), (-3,0)

    std::cout << "Computing ellipse (without randomization)...";
    std::cout.flush();
    Min_ellipse  me1( P, P+n, false);    // very slow
    std::cout << "done." << std::endl;

    std::cout << "Computing ellipse (with randomization)...";
    std::cout.flush();
    Min_ellipse  me2( P, P+n, true);     // fast
    std::cout << "done." << std::endl;

    // because all input points are collinear, the ellipse is
    // degenerate and equals a line segment; the ellipse has
    // two support points
    assert(me2.is_degenerate());
    assert(me2.number_of_support_points()==2);

    // prettyprinting
    CGAL::set_pretty_mode( std::cout);
    std::cout << me2;

    // in general, the ellipse is not explicitly representable
    // over the input number type NT; when you use the default
    // traits class CGAL::Min_ellipse_2_traits_2<K>, you can
    // get double approximations for the coefficients of the
    // underlying conic curve. NOTE: this curve only exists
    // in the nondegenerate case!

    me2.insert(Point(0,1)); // resolves the degeneracy
    assert(!me2.is_degenerate());

    // get the coefficients
    double r,s,t,u,v,w;
    me2.ellipse().double_coefficients( r, s, t, u, v, w);
    std::cout << "ellipse has the equation " <<
      r << " x^2 + " <<
      s << " y^2 + " <<
      t << " xy + " <<
      u << " x + " <<
      v << " y + " <<
      w << " = 0." << std::endl;

    return 0;
}

