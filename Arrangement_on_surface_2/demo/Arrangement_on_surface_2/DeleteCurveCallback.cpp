#include "DeleteCurveCallback.h"
#include "ArrangementTypes.h"

#if 0
template < >
double
DeleteCurveCallback< Seg_arr >::Compute_squared_distance_2::
operator() ( const Point& p, const X_monotone_curve_2& c )
{
    std::cout << "seg_arr curve point distance stub " << std::endl;
    Point p1 = c.source( );
    Point p2 = c.target( );
    Segment seg( p1, p2 );

    return CGAL::to_double( CGAL::squared_distance( p, seg ) );
}
#endif

