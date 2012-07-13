#ifndef CGAL_ARRANGEMENTS_DEMO_UTILS_H
#define CGAL_ARRANGEMENTS_DEMO_UTILS_H
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/iterator.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsView>
#include <QGraphicsScene>
#include "ArrangementTypes.h"
#include <CGAL/Arr_walk_along_line_point_location.h>

class QGraphicsScene;

template < class ArrTraits >
class ArrTraitsAdaptor
{ };

template < class Kernel_ >
class ArrTraitsAdaptor< CGAL::Arr_segment_traits_2< Kernel_ > >
{
public:
    typedef Kernel_ Kernel;
    typedef CGAL::Arr_segment_traits_2< Kernel > ArrTraits;
    typedef typename ArrTraits::Point_2 Point_2;
};

template < class SegmentTraits >
class ArrTraitsAdaptor< CGAL::Arr_polyline_traits_2< SegmentTraits > >
{
public:
    typedef CGAL::Arr_polyline_traits_2< SegmentTraits > ArrTraits;
    typedef typename SegmentTraits::Kernel Kernel;
    typedef typename ArrTraits::Point_2 Point_2;
};

template < class RatKernel, class AlgKernel, class NtTraits >
class ArrTraitsAdaptor< CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > >
{
public:
    typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > ArrTraits;
    typedef AlgKernel Kernel;
    typedef typename ArrTraits::Point_2 Point_2;
};

template < class ArrTraits >
class Compute_squared_distance_2_base
{
public:
    typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel Kernel;
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
    typedef Kernel_ Kernel;
    typedef CGAL::Arr_segment_traits_2< Kernel > Traits;
    typedef Compute_squared_distance_2_base< Traits > Superclass;
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

template < class Kernel_ >
class Compute_squared_distance_2< CGAL::Arr_polyline_traits_2< Kernel_ > > :
    public Compute_squared_distance_2_base< CGAL::Arr_polyline_traits_2< Kernel_ > >
{
public:
    typedef Kernel_ Kernel;
    typedef CGAL::Arr_polyline_traits_2< Kernel > Traits;
    typedef Compute_squared_distance_2_base< Traits > Superclass;
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename Curve_2::const_iterator Curve_const_iterator;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

    FT operator() ( const Point_2& p, const X_monotone_curve_2& c ) const
    {
        Curve_const_iterator ps = c.begin();
        Curve_const_iterator pt = ps; pt++;
        bool first = true;
        FT min_dist = 0;

        while ( pt != c.end() )
        {
            const Point_2& source = *ps;
            const Point_2& target = *pt;
            Segment_2 seg( source, target );
            FT dist = this->squared_distance( p, seg );

            if ( first || dist < min_dist )
            {
                first = false;
                min_dist = dist;
            }
            ps++; pt++;
        }

        return min_dist;
    }
};

template < class RatKernel, class AlgKernel, class NtTraits >
class Compute_squared_distance_2< CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > > :
    public Compute_squared_distance_2_base< CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > >
{
public:
    typedef AlgKernel Kernel;
    typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
    typedef Compute_squared_distance_2_base< Traits > Superclass;
    typedef typename Traits::Point_2 Conic_point_2; // _Conic_point_2< AlgKernel > : public AlgKernel::Point_2
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;

    FT operator() ( const Point_2& p, const X_monotone_curve_2& c ) const
    {
        // Get the co-ordinates of the curve's source and target.
        double sx = CGAL::to_double( c.source( ).x( ) ),
               sy = CGAL::to_double( c.source( ).y( ) ),
               tx = CGAL::to_double( c.target( ).x( ) ),
               ty = CGAL::to_double( c.target( ).y( ) );

        if ( c.orientation( ) == CGAL::COLLINEAR )
        {
            Point_2 ps = c.source( );
            Point_2 pt = c.target( );
            Segment_2 seg( ps, pt );

            return CGAL::squared_distance( p, seg );
        }
        else
        {
            // If the curve is monotone, than its source and its target has the
            // extreme x co-ordinates on this curve.
            bool is_source_left = (sx < tx);
            //int  x_min = is_source_left ? (*w).x_pixel(sx) : (*w).x_pixel(tx);
            //int  x_max = is_source_left ? (*w).x_pixel(tx) : (*w).x_pixel(sx);
            //double   prev_x = is_source_left ? sx : tx;
            //double   prev_y = is_source_left ? sy : ty;
            //double   curr_x, curr_y;
            //int      x;
            //Arr_conic_point_2 px;

            bool first = true;
            FT min_dist( 100000000 );
            AlgKernel ker;

            int n = 100; // TODO: get an adaptive approximation
            std::pair< double, double >* app_pts = new std::pair< double, double >[ n + 1 ];
            std::pair< double, double >* end_pts = c.polyline_approximation( n, app_pts );
            std::pair< double, double >* p_curr = app_pts;
            std::pair< double, double >* p_next = p_curr + 1;
            do
            {
                Point_2 p1( p_curr->first, p_curr->second );
                Point_2 p2( p_next->first, p_next->second );
                Segment_2 seg( p1, p2 );

                FT dist = CGAL::squared_distance( p, seg );
                if ( first || dist < min_dist )
                {
                    first = true;
                    min_dist = dist;
                }

                p_curr++;
                p_next++;
            } while ( p_next != end_pts );

            return min_dist;
        }
    }
};

// TODO: Make Construct_x_monotone_subcurve_2 more generic
template < class ArrTraits >
class Construct_x_monotone_subcurve_2
{
public:
    typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel Kernel;
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

    Construct_x_monotone_subcurve_2( ):
        intersect_2( this->traits.intersect_2_object( ) ),
        split_2( this->traits.split_2_object( ) ),
        compare_x_2( this->traits.compare_x_2_object( ) ),
        construct_x_monotone_curve_2( this->traits.construct_x_monotone_curve_2_object( ) ),
        construct_min_vertex_2( this->traits.construct_min_vertex_2_object( ) ),
        construct_max_vertex_2( this->traits.construct_max_vertex_2_object( ) )
    { }

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
            // FIXME: handle vertical lines properly
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
    ArrTraits traits;
    Intersect_2 intersect_2;
    Split_2 split_2;
    Compare_x_2 compare_x_2;
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
    Construct_min_vertex_2 construct_min_vertex_2;
    Construct_max_vertex_2 construct_max_vertex_2;
}; // class Construct_x_monotone_subcurve_2

template < class RatKernel, class AlgKernel, class NtTraits >
class Construct_x_monotone_subcurve_2< CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > >
{
public:
    typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > ArrTraits;
    typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename AlgKernel::Point_2 Point_2;

    /*
    Return the subcurve of curve bracketed by pLeft and pRight.
    */
    X_monotone_curve_2 operator() ( const X_monotone_curve_2& curve, const Point_2& pLeft, const Point_2& pRight )
    {
        // find the points on the curve
        Point_2 left = curve.point_at_x( pLeft );
        Point_2 right = curve.point_at_x( pRight );

        // make sure the points are oriented in the direction that the curve is going
        AlgKernel ker;
        if (! (((curve.is_directed_right( )) &&
            ker.compare_xy_2_object() ( left, right ) == CGAL::SMALLER) ||
            ((! curve.is_directed_right( )) &&
            ker.compare_xy_2_object() ( left, right ) == CGAL::LARGER)))
        {
            Point_2 tmp = left;
            left = right;
            right = tmp;
        }

        X_monotone_curve_2 res = curve.trim( left, right );
        return res;
    }
}; // class Construct_x_monotone_subcurve_2 for Arr_conic_traits_2


template < class K_ >
class SnapStrategy
{
public:
    typedef K_ Kernel;
    typedef typename Kernel::Point_2 Point_2;

    virtual Point_2 snapPoint( QGraphicsSceneMouseEvent* event ) = 0;
    void setScene( QGraphicsScene* scene_ );

protected:
    SnapStrategy( QGraphicsScene* scene_ );
    QRectF viewportRect( ) const;

    QGraphicsScene* scene;
}; // class SnapStrategy

template < class K_ >
SnapStrategy< K_ >::
SnapStrategy( QGraphicsScene* scene_ ):
    scene( scene_ )
{ }

template < class K_ >
QRectF 
SnapStrategy< K_ >::
viewportRect( ) const
{
    QRectF res;
    if ( this->scene == NULL )
    {
        return res;
    }

    QList< QGraphicsView* > views = this->scene->views( );
    if ( views.size( ) == 0 )
    {
        return res;
    }
    // assumes the first view is the right one
    QGraphicsView* viewport = views.first( );
    QPointF p1 = viewport->mapToScene( 0, 0 );
    QPointF p2 = viewport->mapToScene( viewport->width( ), viewport->height( ) );
    res = QRectF( p1, p2 );

    return res;
}

template < class K_ >
void
SnapStrategy< K_ >::
setScene( QGraphicsScene* scene_ )
{
    this->scene = scene_;
}

template < class K_ >
class SnapToGridStrategy : public SnapStrategy< K_ >
{
public:
    typedef K_ Kernel;
    typedef typename Kernel::Point_2 Point_2;
    typedef SnapStrategy< Kernel > Superclass;

    SnapToGridStrategy( ):
        Superclass( NULL ),
        gridSize( 50 )
    { }

    SnapToGridStrategy( QGraphicsScene* scene ):
        Superclass( scene ),
        gridSize( 50 )
    { }

    Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
    {
        QPointF clickedPoint = event->scenePos( );
        QRectF viewportRect = this->viewportRect( );
        if ( viewportRect == QRectF( ) )
        {
            return this->convert( event->scenePos( ) );
        }

        qreal d( this->gridSize / 2.0 );
        int left = int( viewportRect.left( ) ) - (int( viewportRect.left( ) ) % this->gridSize);
        int right = int( viewportRect.right( ) ) + (this->gridSize - int( viewportRect.right( ) ) % this->gridSize);
        int x = clickedPoint.x( );
        int y = clickedPoint.y( );
        for ( int i = left - this->gridSize; i <= right; i += this->gridSize )
        {
            if ( i - d <= clickedPoint.x( ) && clickedPoint.x( ) <= i + d )
            {
                x = i;
                break;
            }
        }
        int top = int( viewportRect.top( ) ) - (int( viewportRect.top( ) ) % this->gridSize);
        int bottom = int( viewportRect.bottom( ) ) + (this->gridSize - int( viewportRect.bottom( ) ) % this->gridSize);
        for ( int i = top - this->gridSize; i <= bottom; i += this->gridSize )
        {
            if ( i - d <= clickedPoint.y( ) && clickedPoint.y( ) <= i + d )
            {
                y = i;
                break;
            }
        }
        return this->convert( QPointF( x, y ) );
    }

    void setGridSize( int size )
    {
        this->gridSize = size;
    }

protected:
    int gridSize;
    CGAL::Qt::Converter< Kernel > convert;
}; // class SnapToGridStrategy

template < class Arr_ >
class SnapToArrangementVertexStrategy:
    public SnapStrategy< typename ArrTraitsAdaptor< typename Arr_::Geometry_traits_2 >::Kernel >
{
public:
    typedef Arr_ Arrangement;
    typedef typename Arrangement::Geometry_traits_2 Traits;
    typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
    typedef SnapStrategy< Kernel > Superclass;
    typedef typename Arrangement::Vertex_iterator Vertex_iterator;
    typedef typename Kernel::Compute_squared_distance_2 Compute_squared_distance_2;
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_2 Point_2;

    SnapToArrangementVertexStrategy( ):
        Superclass( NULL ),
        arrangement( NULL )
    { }

    SnapToArrangementVertexStrategy( Arrangement* arr, QGraphicsScene* scene_ ):
        Superclass( scene_ ),
        arrangement( arr )
    { }

    Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
    {
        Point_2 clickedPoint = this->convert( event->scenePos( ) );
        Point_2 closestPoint = clickedPoint;
        bool first = true;
        FT minDist( 0 );
        QRectF viewportRect = this->viewportRect( );
        if ( viewportRect == QRectF( ) )
        {
            return this->convert( event->scenePos( ) );
        }

        FT maxDist( ( viewportRect.right( ) - viewportRect.left( ) ) / 4.0 );
        for ( Vertex_iterator vit = this->arrangement->vertices_begin( ); 
                vit != this->arrangement->vertices_end( ); ++vit )
        {
            Point_2 point = vit->point( );
            FT dist = this->compute_squared_distance_2( clickedPoint, point );
            if ( first || ( dist < minDist ) )
            {
                first = false;
                minDist = dist;
                closestPoint = point;
            }
        }
        if ( ! first && minDist < maxDist )
        {
            return closestPoint;
        }
        else
        {
            return this->convert( event->scenePos( ) );
        }
    }

    void setArrangement( Arrangement* arr )
    {
        this->arrangement = arr;
    }

protected:
    Arrangement* arrangement;
    Compute_squared_distance_2 compute_squared_distance_2;
    CGAL::Qt::Converter< Kernel > convert;
}; // class SnapToArrangementVertexStrategy

template < class Arr_ >
class Find_nearest_edge
{
public: // typedefs
    typedef Arr_ Arrangement;
    typedef typename Arrangement::Geometry_traits_2 ArrTraits;
    typedef Compute_squared_distance_2< ArrTraits > Point_curve_distance;
    typedef typename ArrTraits::X_monotone_curve_2 X_monotone_curve_2;
    typedef CGAL::Arr_walk_along_line_point_location< Arrangement > Point_location_strategy;
    typedef typename ArrTraitsAdaptor< ArrTraits >::Kernel Kernel;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Arrangement::Face_const_handle Face_const_handle;
    typedef typename Arrangement::Halfedge_const_handle Halfedge_const_handle;
    typedef typename Arrangement::Vertex_const_handle Vertex_const_handle;
    typedef typename Arrangement::Ccb_halfedge_const_circulator Ccb_halfedge_const_circulator;
    typedef typename Point_curve_distance::FT FT;
    typedef typename Arrangement::Hole_const_iterator Hole_const_iterator;
    typedef typename Arrangement::Halfedge_around_vertex_const_circulator Halfedge_around_vertex_const_circulator;

public: // constructors
    Find_nearest_edge( Arrangement* arr_ ):
        arr( arr_ ),
        pointLocationStrategy( Point_location_strategy( *arr_ ) )
    { }

public: // member methods
    Halfedge_const_handle operator()( const Point_2& queryPt )
    {
        CGAL::Object pointLocationResult = this->pointLocationStrategy.locate( queryPt );
        Face_const_handle face = this->getFace( pointLocationResult );
        bool first = 1;
        X_monotone_curve_2 closestCurve;
        Halfedge_const_handle closestEdge;
        FT minDist( 0 );

        if ( ! face->is_unbounded( ) )
        { // it is an interior face so it has a ccb
            Ccb_halfedge_const_circulator cc = face->outer_ccb( );
            do
            {
                X_monotone_curve_2 curve = cc->curve( );
                FT dist = this->pointCurveDistance( queryPt, curve );
                if ( first || dist < minDist )
                {
                    first = 0;
                    minDist = dist;
                    closestEdge = cc;
                }
            }
            while ( ++cc != face->outer_ccb( ) );
        }
        Hole_const_iterator hit; 
        Hole_const_iterator eit = face->holes_end( );
        for ( hit = face->holes_begin( ); hit != eit; ++hit )
        { // check any holes inside this face
            Ccb_halfedge_const_circulator cc = *hit;
            do
            {
                X_monotone_curve_2 curve = cc->curve( );
                FT dist = this->pointCurveDistance( queryPt, curve );
                if ( first || dist < minDist )
                {
                    first = 0;
                    minDist = dist;
                    closestEdge = cc;
                }
                cc++;
            }
            while ( cc != *hit );
        }

        return closestEdge;
    }

protected: // member methods
    Face_const_handle getFace( const CGAL::Object& obj )
    {
        Face_const_handle f;
        if ( CGAL::assign( f, obj ) )
            return f;

        Halfedge_const_handle he;
        if (CGAL::assign( he, obj ))
            return (he->face( ));

        Vertex_const_handle v;
        CGAL_assertion(CGAL::assign( v, obj ));
        CGAL::assign( v, obj );
        if ( v->is_isolated( ) )
            return v->face( );
        Halfedge_around_vertex_const_circulator eit = v->incident_halfedges( );
        return  (eit->face( ));
    }

protected: // member fields
    Arrangement* arr;
    Point_curve_distance pointCurveDistance;
    Point_location_strategy pointLocationStrategy;

}; // class Find_nearest_edge

#endif // CGAL_ARRANGEMENTS_DEMO_UTILS_H
