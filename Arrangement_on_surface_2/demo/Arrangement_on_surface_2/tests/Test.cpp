#include <CGAL/Exact_circular_kernel_2.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_circular_arc_traits_2.h>

#include <fstream>

typedef CGAL::Exact_circular_kernel_2 CircularKernel;
typedef CGAL::Arr_circular_arc_traits_2< CircularKernel > Traits;
typedef CGAL::Arr_default_dcel< Traits > Dcel;
typedef CGAL::Arrangement_with_history_2< Traits, Dcel > Arrangement;

int main( )
{
    Arrangement arr;

    std::ifstream ifs( "curves" );
    Traits::Curve_2 testArc, testArc2;
    ifs >> testArc >> testArc2;
    ifs.close( );

    std::cout << testArc << std::endl;
    std::cout << testArc2 << std::endl;

    CGAL::insert( arr, testArc );
    CGAL::insert( arr, testArc2 );

    return 0;
}
