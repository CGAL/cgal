#ifndef CGAL_ARRANGEMENTS_DEMO_UTILS_H
#define CGAL_ARRANGEMENTS_DEMO_UTILS_H
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/iterator.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsView>
#include <QGraphicsScene>

class QGraphicsScene;

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

    Intersect_2 intersect_2;
    Split_2 split_2;
    Compare_x_2 compare_x_2;
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
    Construct_min_vertex_2 construct_min_vertex_2;
    Construct_max_vertex_2 construct_max_vertex_2;
};

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
    public SnapStrategy< typename Arr_::Geometry_traits_2::Kernel >
{
public:
    typedef Arr_ Arrangement;
    typedef typename Arrangement::Geometry_traits_2 Traits;
    typedef typename Traits::Kernel Kernel;
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

#endif // CGAL_ARRANGEMENTS_DEMO_UTILS_H
