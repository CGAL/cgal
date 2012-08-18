#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_circular_arc_traits_2.h>

#include <fstream>

typedef CGAL::Exact_circular_kernel_2 CircularKernel;
typedef CGAL::Arr_circular_arc_traits_2< CircularKernel > Traits;
typedef CGAL::Arr_default_dcel< Traits > Dcel;
typedef CGAL::Arrangement_with_history_2< Traits, Dcel > Arrangement;
typedef Traits::X_monotone_curve_2 Circular_arc_2;
typedef CircularKernel::Circle_2 Circle_2;
typedef CircularKernel::Point_2 Point_2;

typedef CGAL::Arrangement_2< Traits, Dcel > Arrangement_without_history;

int main( )
{
    Arrangement arr;
    Arrangement_without_history arr2;

    std::ifstream ifs( "curves" );
    Traits::Curve_2 testArc, testArc2;
    ifs >> testArc >> testArc2;
    ifs.close( );

    std::cout << testArc << std::endl;
    std::cout << testArc2 << std::endl;

    CGAL::insert( arr, testArc );
    CGAL::insert( arr, testArc2 );

    Point_2 p1( 0, 0 );
    Point_2 p2( 1, 0 );
    Point_2 p3( 0, 1 );
    Circle_2 circle( p1, p2, p3 );
    Circular_arc_2 arc( circle );
    Traits::Curve_2 curveBox = arc;
    
    std::vector< Circular_arc_2 > arcs;
    arcs.push_back( arc );
    CGAL::insert( arr2, arcs.begin( ), arcs.end( ) );

    return 0;
}
