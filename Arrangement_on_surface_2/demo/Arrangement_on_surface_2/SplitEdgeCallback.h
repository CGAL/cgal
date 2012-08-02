#ifndef SPLIT_EDGE_CALLBACK_H
#define SPLIT_EDGE_CALLBACK_H
#include "Callback.h"
#include <QEvent>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <CGAL/Qt/Converter.h>
#include "CurveGraphicsItem.h"
#include <CGAL/Arrangement_with_history_2.h>
#include "Utils.h"
#include "ISnappable.h"

class SplitEdgeCallbackBase : public CGAL::Qt::Callback, public ISnappable
{
public:
    void setSnappingEnabled( bool b );
    void setSnapToGridEnabled( bool b );

protected:
    SplitEdgeCallbackBase( QObject* parent );

    bool snappingEnabled;
    bool snapToGridEnabled;
}; // SplitEdgeCallbackBase

/**
Handles splitting of arrangement curves selected from the scene.

The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < class Arr_ >
class SplitEdgeCallback : public SplitEdgeCallbackBase
{
public:
    typedef Arr_ Arrangement;
    typedef typename Arrangement::Halfedge_handle Halfedge_handle;
    typedef typename Arrangement::Halfedge_iterator Halfedge_iterator;
    typedef typename Arrangement::Vertex_iterator Vertex_iterator;
    typedef typename Arrangement::Geometry_traits_2 Traits;
    typedef typename Arrangement::Curve_handle Curve_handle;
    typedef typename Arrangement::Originating_curve_iterator Originating_curve_iterator;
    typedef typename Arrangement::Induced_edge_iterator Induced_edge_iterator;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
    typedef typename Traits::Intersect_2 Intersect_2;
    typedef typename Traits::Equal_2 Equal_2;
    typedef typename Traits::Multiplicity Multiplicity;
    typedef typename ArrTraitsAdaptor< Traits >::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Kernel::FT FT;

    SplitEdgeCallback( Arrangement* arr_, QObject* parent );
    void setScene( QGraphicsScene* scene_ );
    void reset( );
    
    void slotModelChanged( );

protected:
    void mousePressEvent( QGraphicsSceneMouseEvent *event );
    void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
    virtual Point_2 snapPoint( QGraphicsSceneMouseEvent *event );

    template < class TTraits >
    void splitEdges( const Point_2& pt, TTraits traits );

    template < class CircularKernel >
    void splitEdges( const Point_2& pt, CGAL::Arr_circular_arc_traits_2< CircularKernel > traits );

    template < class TTraits >
    void updateGuide( const Point_2& pt, TTraits traits );

    template < class CircularKernel >
    void updateGuide( const Point_2& pt, CGAL::Arr_circular_arc_traits_2< CircularKernel > traits );

    Traits traits;
    CGAL::Qt::Converter< Kernel > convert;
    Arrangement* arr;
    bool hasFirstPoint;
    Point_2 p1;
    Point_2 p2;
    QGraphicsLineItem segmentGuide;

    Intersect_2 intersectCurves;
    Equal_2 areEqual;
    SnapToArrangementVertexStrategy< Arrangement > snapToVertexStrategy;
    SnapToGridStrategy< Kernel > snapToGridStrategy;
}; // class SplitEdgeCallback

template < class Arr_ >
SplitEdgeCallback< Arr_ >::
SplitEdgeCallback( Arrangement* arr_, QObject* parent ):
    SplitEdgeCallbackBase( parent ),
    arr( arr_ ),
    hasFirstPoint( false ),
    intersectCurves( this->traits.intersect_2_object( ) ),
    areEqual( this->traits.equal_2_object( ) )
{
    this->snapToVertexStrategy.setArrangement( arr_ );

    QObject::connect( this, SIGNAL( modelChanged( ) ),
        this, SLOT( slotModelChanged( ) ) );
}

template < class Arr_ >
void 
SplitEdgeCallback< Arr_ >::
setScene( QGraphicsScene* scene_ )
{
    this->scene = scene_;
    this->snapToVertexStrategy.setScene( scene_ );
    this->snapToGridStrategy.setScene( scene_ );
    if ( this->scene )
    {
        this->scene->addItem( &( this->segmentGuide ) );
    }
}

template < class Arr_ >
void
SplitEdgeCallback< Arr_ >::
reset( )
{
    this->hasFirstPoint = false;
    this->segmentGuide.setLine( 0, 0, 0, 0 );
    emit modelChanged( );
}

template < class Arr_ >
void
SplitEdgeCallback< Arr_ >::
slotModelChanged( )
{
    this->segmentGuide.update( );
}

template < class Arr_ >
void 
SplitEdgeCallback< Arr_ >::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{
    Point_2 clickedPoint = this->snapPoint( event );
    this->splitEdges( clickedPoint, Traits( ) );
}

template < class Arr_ >
template < class TTraits >
void
SplitEdgeCallback< Arr_ >::
splitEdges( const Point_2& clickedPoint, TTraits traits )
{
    typename TTraits::Construct_x_monotone_curve_2 construct_x_monotone_curve_2 =
        traits.construct_x_monotone_curve_2_object( );
    if ( ! this->hasFirstPoint )
    {
        this->p1 = clickedPoint;
        this->hasFirstPoint = true;
    }
    else
    {
        this->p2 = clickedPoint;
        X_monotone_curve_2 splitCurve =
            construct_x_monotone_curve_2( this->p1, this->p2 );
        for ( Halfedge_iterator hei = this->arr->halfedges_begin( );
            hei != this->arr->halfedges_end( ); ++hei )
        {
            X_monotone_curve_2 curve = hei->curve( );
            CGAL::Object res;
            CGAL::Oneset_iterator< CGAL::Object > oi( res );
            this->intersectCurves( splitCurve, curve, oi );
            std::pair< Point_2, Multiplicity > pair;
            if ( hei == this->arr->halfedges_end( ) )
                continue;
            if ( CGAL::assign( pair, res ) )
            {
                Point_2 splitPoint = pair.first;
                if ( ( ! hei->source( )->is_at_open_boundary( ) && this->areEqual( hei->source( )->point( ), splitPoint ) ) ||
                    ( ! hei->target( )->is_at_open_boundary( ) && this->areEqual( hei->target( )->point( ), splitPoint ) ) )
                {
                    continue;
                }
                this->arr->split_edge( hei, splitPoint );
            }
        }

        this->reset( );
    }

    emit modelChanged( );
}

template < class Arr_ >
template < class CircularKernel >
void
SplitEdgeCallback< Arr_ >::
splitEdges( const Point_2& clickedPoint, CGAL::Arr_circular_arc_traits_2< CircularKernel > traits )
{
    std::cout << "Circular arc split edges stub" << std::endl;
}

template < class Arr_ >
void 
SplitEdgeCallback< Arr_ >::
mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{
    Point_2 clickedPoint = this->snapPoint( event );
    this->updateGuide( clickedPoint, this->traits );
#if 0
    if ( this->hasFirstPoint )
    { // provide visual feedback for where the split line is
        Point_2 currentPoint = clickedPoint;
        Segment_2 currentSegment( this->p1, currentPoint );
        QLineF qSegment = this->convert( currentSegment );
        this->segmentGuide.setLine( qSegment );
        emit modelChanged( );
    }
#endif
}

template < class Arr_ >
template < class TTraits >
void
SplitEdgeCallback< Arr_ >::
updateGuide( const Point_2& clickedPoint, TTraits traits )
{
    if ( this->hasFirstPoint )
    { // provide visual feedback for where the split line is
        Point_2 currentPoint = clickedPoint;
        Segment_2 currentSegment( this->p1, currentPoint );
        QLineF qSegment = this->convert( currentSegment );
        this->segmentGuide.setLine( qSegment );
        emit modelChanged( );
    }
}

template < class Arr_ >
template < class CircularKernel >
void
SplitEdgeCallback< Arr_ >::
updateGuide( const Point_2& clickedPoint, CGAL::Arr_circular_arc_traits_2< CircularKernel > traits )
{
    if ( this->hasFirstPoint )
    { // provide visual feedback for where the split line is
        Point_2 currentPoint = clickedPoint;
        typename CircularKernel::Point_2 pt1( CGAL::to_double( this->p1.x( ) ), CGAL::to_double( this->p1.y( ) ) );
        typename CircularKernel::Point_2 pt2( CGAL::to_double( currentPoint.x( ) ), CGAL::to_double( currentPoint.y( ) ) );
        Segment_2 currentSegment( pt1, pt2 );
        QLineF qSegment = this->convert( currentSegment );
        this->segmentGuide.setLine( qSegment );
        emit modelChanged( );
    }
}

template < class Arr_ >
typename SplitEdgeCallback< Arr_ >::Point_2
SplitEdgeCallback< Arr_ >::
snapPoint( QGraphicsSceneMouseEvent* event )
{
    if ( this->snapToGridEnabled )
    {
        return this->snapToGridStrategy.snapPoint( event );
    }
    if ( this->snappingEnabled )
    {
        return this->snapToVertexStrategy.snapPoint( event );
    }
    else
    {
        return this->convert( event->scenePos( ) );
    }
}

#endif // SPLIT_EDGE_CALLBACK_H
