#ifndef CGAL_QT_GRAPHICS_VIEW_CURVE_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_CURVE_INPUT_H
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QEvent>
#include <QGraphicsLineItem>
#include <QGraphicsSceneMouseEvent>
#include <iostream>
#include "Callback.h"
#include "ISnappable.h"

namespace CGAL {
namespace Qt {

class GraphicsViewCurveInputBase:
    public GraphicsViewInput, public ISnappable
{
public:
    virtual void setScene( QGraphicsScene* scene_ );
    QGraphicsScene* getScene( ) const;

    void setSnappingEnabled( bool b );
    void setSnapToGridEnabled( bool b );

protected:
    GraphicsViewCurveInputBase( QObject* parent );
    virtual void mouseMoveEvent( QGraphicsSceneMouseEvent* event );
    virtual void mousePressEvent( QGraphicsSceneMouseEvent* event );
    virtual bool eventFilter( QObject* obj, QEvent* event );

    QRectF viewportRect( ) const;

    QGraphicsScene* scene;
    bool snappingEnabled;
    bool snapToGridEnabled;

}; // class GraphicsViewCurveInputBase

template < class ArrTraits >
class GraphicsViewCurveInput:
    public GraphicsViewCurveInputBase
{ };

/**
Specialization of GraphicsViewCurveInput for Arr_segment_traits_2; handles
user-guided generation of line segment curves.
*/
template < class Kernel_ >
class GraphicsViewCurveInput< CGAL::Arr_segment_traits_2< Kernel_ > >:
    public GraphicsViewCurveInputBase
{
public:
    typedef Kernel_ Kernel;
    typedef CGAL::Arr_segment_traits_2< Kernel > Traits;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;

    GraphicsViewCurveInput( QObject* parent ):
        GraphicsViewCurveInputBase( parent ),
        second( false )
    { }

protected:
    void mouseMoveEvent( QGraphicsSceneMouseEvent* event )
    {
        if ( this->second )
        {
            Point_2 clickedPoint = this->snapPoint( event );
            Segment_2 segment( this->p1, clickedPoint );
            QLineF qSegment = this->convert( segment );
            this->segmentGuide.setLine( qSegment );
        }
    }

    void mousePressEvent( QGraphicsSceneMouseEvent* event )
    {
        if ( !this->second )
        {
            this->second = true;
            this->p1 = this->snapPoint( event );
            QPointF pt = this->convert( this->p1 );
            this->segmentGuide.setLine( pt.x( ), pt.y( ), pt.x( ), pt.y( ) );
            if ( this->scene != NULL )
            {
                this->scene->addItem( &( this->segmentGuide ) );
            }
        }
        else
        {
            this->second = false;
            this->p2 = this->snapPoint( event );
            if ( this->scene != NULL )
            {
                this->scene->removeItem( &( this->segmentGuide ) );
            }
            Curve_2 res( this->p1, this->p2 );
            emit generate( CGAL::make_object( res ) );
        }
    }

    // override this to snap to the points you like
    virtual Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
    {
        Point_2 clickedPoint = this->convert( event->scenePos( ) );
        return clickedPoint;
    }

    Converter< Kernel > convert;
    Point_2 p1;
    Point_2 p2;
    bool second;

    QGraphicsLineItem segmentGuide;
}; // class GraphicsViewCurveInput< CGAL::Arr_segment_traits_2< Kernel_ > >

/**
Specialization of GraphicsViewCurveInput for Arr_polyline_traits_2; handles
user-guided generation of line segment curves.
*/
template < class SegmentTraits >
class GraphicsViewCurveInput< CGAL::Arr_polyline_traits_2< SegmentTraits > >:
    public GraphicsViewCurveInputBase
{
public:
    typedef CGAL::Arr_polyline_traits_2< SegmentTraits > Traits;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename SegmentTraits::Kernel Kernel;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;

    GraphicsViewCurveInput( QObject* parent ):
        GraphicsViewCurveInputBase( parent )
    { }

protected:
    void mouseMoveEvent( QGraphicsSceneMouseEvent* event )
    {
        if ( ! this->polylineGuide.empty( ) )
        {
            Point_2 clickedPoint = this->snapPoint( event );
            // TODO: make it work for the latest line segment
            Segment_2 segment( this->points.back( ), clickedPoint );
            QLineF qSegment = this->convert( segment );
            this->polylineGuide.back( )->setLine( qSegment );
        }
    }

    void mousePressEvent( QGraphicsSceneMouseEvent* event )
    {
        Point_2 clickedPoint = this->snapPoint( event );
        if ( this->points.empty( ) )
        { // first
            // add clicked point to polyline
            this->points.push_back( clickedPoint );

            QPointF pt = this->convert( clickedPoint );
            QGraphicsLineItem* lineItem = new QGraphicsLineItem( pt.x( ), pt.y( ), pt.x( ), pt.y( ) );
            this->polylineGuide.push_back( lineItem );
            if ( this->scene != NULL )
            {
                this->scene->addItem( this->polylineGuide.back( ) );
            }
        }
        else
        {
            // add clicked point to polyline
            this->points.push_back( clickedPoint );

            if ( event->button( ) == ::Qt::RightButton )
            { // finalize polyline input
                for ( int i = 0; i < this->polylineGuide.size( ); ++i )
                {
                    if ( this->scene != NULL )
                    {
                        this->scene->removeItem( this->polylineGuide[ i ] );
                    }
                    delete this->polylineGuide[ i ];
                }
                this->polylineGuide.clear( );
                Curve_2 res( this->points.begin( ), this->points.end( ) );
                this->points.clear( );

                emit generate( CGAL::make_object( res ) );
            }
            else
            { // start the next segment
                QPointF pt = this->convert( clickedPoint );
                QGraphicsLineItem* lineItem = new QGraphicsLineItem( pt.x( ), pt.y( ), pt.x( ), pt.y( ) );
                this->polylineGuide.push_back( lineItem );
                if ( this->scene != NULL )
                {
                    this->scene->addItem( this->polylineGuide.back( ) );
                }
            }
        }
    }

    // override this to snap to the points you like
    virtual Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
    {
        Point_2 clickedPoint = this->convert( event->scenePos( ) );
        return clickedPoint;
    }

    Converter< Kernel > convert;
    std::vector< Point_2 > points;

    std::vector< QGraphicsLineItem* > polylineGuide;
}; // class GraphicsViewCurveInput< CGAL::Arr_polyline_traits_2< SegmentTraits > >

/**
Specialization of GraphicsViewCurveInput for Arr_conic_traits_2; handles
user-guided generation of conic curves.
*/
template < class RatKernel, class AlgKernel, class NtTraits >
class GraphicsViewCurveInput< CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > >:
    public GraphicsViewCurveInputBase
{
public:
    typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
    typedef typename Traits::Curve_2 Curve_2;
    typedef typename Traits::Point_2 Point_2; // basically, an AlgKernel::Point_2 with metadata
    typedef typename Traits::Construct_x_monotone_curve_2 Construct_x_monotone_curve_2;
    typedef AlgKernel Kernel;
    //typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename RatKernel::Point_2 Rat_point_2;
    typedef typename RatKernel::Segment_2 Rat_segment_2;

    GraphicsViewCurveInput( QObject* parent ):
        GraphicsViewCurveInputBase( parent ),
        construct_x_monotone_curve_2( this->traits.construct_x_monotone_curve_2_object( ) )
    { }

protected:
    void mouseMoveEvent( QGraphicsSceneMouseEvent* event )
    {
        if ( ! this->polylineGuide.empty( ) )
        {
            Point_2 clickedPoint = this->snapPoint( event );
            // TODO: make it work for the latest line segment
            Segment_2 segment( this->points.back( ), clickedPoint );
            QLineF qSegment = this->convert( segment );
            this->polylineGuide.back( )->setLine( qSegment );
        }
    }

    void mousePressEvent( QGraphicsSceneMouseEvent* event )
    {
        Point_2 clickedPoint = this->snapPoint( event );
        this->points.push_back( clickedPoint );

        if ( this->points.size( ) == 1 )
        { // first
            // add clicked point to polyline

            QPointF pt = this->convert( clickedPoint );
            QGraphicsLineItem* lineItem = new QGraphicsLineItem( pt.x( ), pt.y( ), pt.x( ), pt.y( ) );
            this->polylineGuide.push_back( lineItem );
            if ( this->scene != NULL )
            {
                this->scene->addItem( this->polylineGuide.back( ) );
            }
        }
        else
        {
            for ( int i = 0; i < this->polylineGuide.size( ); ++i )
            {
                if ( this->scene != NULL )
                {
                    this->scene->removeItem( this->polylineGuide[ i ] );
                }
                delete this->polylineGuide[ i ];
            }
            this->polylineGuide.clear( );

            //Curve_2 res = this->construct_x_monotone_curve_2( this->points[ 0 ], this->points[ 1 ] );
            double x1 = CGAL::to_double( this->points[ 0 ].x( ) ); 
            double y1 = CGAL::to_double( this->points[ 0 ].y( ) ); 
            double x2 = CGAL::to_double( this->points[ 1 ].x( ) ); 
            double y2 = CGAL::to_double( this->points[ 1 ].y( ) ); 
            Curve_2 res = Curve_2( Rat_segment_2( Rat_point_2( x1, y1 ), Rat_point_2( x2, y2 ) ) );
            std::cout << "res is " << ( (res.is_valid( ))? "" : "not ") << "valid" << std::endl;
            this->points.clear( );

            emit generate( CGAL::make_object( res ) );
        }
    }

    // override this to snap to the points you like
    virtual Point_2 snapPoint( QGraphicsSceneMouseEvent* event )
    {
        Point_2 clickedPoint = this->convert( event->scenePos( ) );
        return clickedPoint;
    }

    Converter< Kernel > convert;
    std::vector< Point_2 > points;
    std::vector< QGraphicsLineItem* > polylineGuide;
    Traits traits;
    Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
}; // class GraphicsViewCurveInput< CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > >

} // namespace Qt
} // namespace CGAL
#endif // CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
