#ifndef CGAL_ARRANGEMENTS_DEMO_UTILS_H
#define CGAL_ARRANGEMENTS_DEMO_UTILS_H
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/iterator.h>

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

template < class ArrTraits >
class Construct_x_monotone_subcurve_2
{
public:
    typedef typename ArrTraits::Kernel Kernel;
    typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename ArrTraits::Split_2 Split_2;
    typedef typename ArrTraits::Intersect_2 Intersect_2;
    typedef typename ArrTraits::Multiplicity Multiplicity;
    typedef typename ArrTraits::Construct_x_monotone_curve_2 Construct_x_monotone_curve_2;
    typedef typename ArrTraits::Construct_min_vertex_2 Construct_min_vertex_2;
    typedef typename ArrTraits::Construct_max_vertex_2 Construct_max_vertex_2;
    typedef typename ArrTraits::Compare_x_2 Compare_x_2;
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Line_2 Line_2;
    typedef typename Kernel::Compute_y_at_x_2 Compute_y_at_x_2;

    /*
    Return the subcurve of curve bracketed by pLeft and pRight.

    We assume pLeft and pRight don't lie on the curve and always do a vertical
    projection.
    */
    X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve, const Point_2& pLeft, const Point_2& pRight )
    {
        Point_2 pMin = this->construct_min_vertex_2( curve );
        Point_2 pMax = this->construct_max_vertex_2( curve );
        X_monotone_curve_2 subcurve;
        X_monotone_curve_2 unusedTrimmings;
        X_monotone_curve_2 finalSubcurve;
        if ( this->compare_x_2( pLeft, pMin ) == CGAL::LARGER )
        {
            CGAL::Bbox_2 c_bbox = curve.bbox( );
            FT splitLineYMin( c_bbox.ymin( ) - 1.0 );
            FT splitLineYMax( c_bbox.ymax( ) + 1.0 );
            Point_2 splitLinePBottom( pLeft.x( ), splitLineYMin );
            Point_2 splitLinePTop( pLeft.x( ), splitLineYMax );
            X_monotone_curve_2 splitLine = 
                this->construct_x_monotone_curve_2( splitLinePBottom, splitLinePTop );
            CGAL::Object res;
            CGAL::Oneset_iterator< CGAL::Object > oi( res );
            this->intersect_2( splitLine, curve, oi );
            std::pair< Point_2, Multiplicity > pair;
            if ( CGAL::assign( pair, res ) )
            {
                Point_2 splitPoint = pair.first;
                this->split_2( curve, splitPoint, unusedTrimmings, subcurve );
            }
        }
        else
        {
            subcurve = curve;
        }
        
        if ( this->compare_x_2( pRight, pMax ) == CGAL::SMALLER )
        {
            CGAL::Bbox_2 c_bbox = subcurve.bbox( );
            FT splitLineYMin( c_bbox.ymin( ) - 1.0 );
            FT splitLineYMax( c_bbox.ymax( ) + 1.0 );
            Point_2 splitLinePBottom( pRight.x( ), splitLineYMin );
            Point_2 splitLinePTop( pRight.x( ), splitLineYMax );
            X_monotone_curve_2 splitLine =
                this->construct_x_monotone_curve_2( splitLinePBottom, splitLinePTop );
            CGAL::Object res;
            CGAL::Oneset_iterator< CGAL::Object > oi( res );
            this->intersect_2( splitLine, subcurve, oi );
            std::pair< Point_2, Multiplicity > pair;
            if ( CGAL::assign( pair, res ) )
            {
                Point_2 splitPoint = pair.first;
                //return X_monotone_curve_2( splitLinePBottom, splitPoint );
                this->split_2( subcurve, splitPoint, finalSubcurve, unusedTrimmings );
            }
        }
        else
        {
            finalSubcurve = subcurve;
        }

        return finalSubcurve;
    }

protected:

    Intersect_2 intersect_2;
    Split_2 split_2;
    Compare_x_2 compare_x_2;
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
    Construct_min_vertex_2 construct_min_vertex_2;
    Construct_max_vertex_2 construct_max_vertex_2;
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
