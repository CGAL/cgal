#ifndef CGAL_ARRANGEMENTS_DEMO_UTILS_H
#define CGAL_ARRANGEMENTS_DEMO_UTILS_H
#include <CGAL/Arr_segment_traits_2.h>

template < class ArrTraits >
class Compute_squared_distance_2_base
{
public:
    typedef typename ArrTraits::Kernel Kernel;
    typedef typename Kernel::FT FT;

    template < class T1, class T2 >
    FT
    operator() ( const T1& t1, const T2& t2 ) const
    {
        return this->squared_distance( t1, t2 );
    }

protected:
    typename Kernel::Compute_squared_distance_2 squared_distance;
};

template < class ArrTraits >
class Compute_squared_distance_2 : public Compute_squared_distance_2_base< ArrTraits >
{ };

template < class Kernel_ >
class Compute_squared_distance_2< CGAL::Arr_segment_traits_2< Kernel_ > > :
    public Compute_squared_distance_2_base< CGAL::Arr_segment_traits_2< Kernel_ > >
{
public:
    typedef CGAL::Arr_segment_traits_2< Kernel_ > Traits;
    typedef Compute_squared_distance_2_base< Traits > Superclass;
    typedef Kernel_ Kernel;
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

    FT operator() ( const Point_2& p, const X_monotone_curve_2& c ) const
    {
        Point_2 p1 = c.source( );
        Point_2 p2 = c.target( );
        Segment_2 seg( p1, p2 );

        return this->squared_distance( p, seg );
    }
};

#if 0
DeleteCurveCallback< TArr >::Compute_squared_distance_2< CGAL::Arr_segment_traits_2< typename TArr::Geometry_traits_2::Kernel > >::
operator() ( const typename TArr::Geometry_traits_2::Kernel::Point& p, const typename TArr::Geometry_traits_2::X_monotone_curve_2& c )
{
    std::cout << "seg_arr curve point distance stub " << std::endl;
    Point p1 = c.source( );
    Point p2 = c.target( );
    Segment seg( p1, p2 );

    return CGAL::to_double( CGAL::squared_distance( p, seg ) );
}
#endif
#endif // CGAL_ARRANGEMENTS_DEMO_UTILS_H
