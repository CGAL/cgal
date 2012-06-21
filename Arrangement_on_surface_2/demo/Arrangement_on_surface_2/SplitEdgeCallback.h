#ifndef SPLIT_EDGE_CALLBACK_H
#define SPLIT_EDGE_CALLBACK_H
#include "Callback.h"
#include <QEvent>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/CurveGraphicsItem.h>
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
template < class TArr >
class SplitEdgeCallback : public SplitEdgeCallbackBase
{
public:
    typedef typename TArr::Halfedge_handle Halfedge_handle;
    typedef typename TArr::Halfedge_iterator Halfedge_iterator;
    typedef typename TArr::Vertex_iterator Vertex_iterator;
    typedef typename TArr::Geometry_traits_2 Traits;
    typedef typename TArr::Curve_handle Curve_handle;
    typedef typename TArr::Originating_curve_iterator Originating_curve_iterator;
    typedef typename TArr::Induced_edge_iterator Induced_edge_iterator;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename Traits::Construct_x_monotone_curve_2 Construct_x_monotone_curve_2;
    typedef typename Traits::Kernel Kernel;
    typedef typename Traits::Intersect_2 Intersect_2;
    typedef typename Traits::Equal_2 Equal_2;
    typedef typename Traits::Multiplicity Multiplicity;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Kernel::FT FT;

    SplitEdgeCallback( TArr* arr_, QObject* parent );
    void setScene( QGraphicsScene* scene_ );
    void reset( );
    
    void slotModelChanged( );

protected:
    void mousePressEvent( QGraphicsSceneMouseEvent *event );
    void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
    virtual Point_2 snapPoint( QGraphicsSceneMouseEvent *event );

    CGAL::Qt::Converter< Kernel > convert;
    QGraphicsScene* scene;
    TArr* arr;
    bool hasFirstPoint;
    Point_2 p1;
    Point_2 p2;
    QGraphicsLineItem segmentGuide;

    Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
    Intersect_2 intersectCurves;
    Equal_2 areEqual;
    SnapToArrangementVertexStrategy< TArr > snapToVertexStrategy;
    SnapToGridStrategy< Kernel > snapToGridStrategy;
}; // class SplitEdgeCallback

template < class TArr >
SplitEdgeCallback< TArr >::
SplitEdgeCallback( TArr* arr_, QObject* parent ):
    SplitEdgeCallbackBase( parent ),
    arr( arr_ ),
    scene( NULL ),
    hasFirstPoint( false )
{
    this->snapToVertexStrategy.setArrangement( arr_ );

    QObject::connect( this, SIGNAL( modelChanged( ) ),
        this, SLOT( slotModelChanged( ) ) );
}

template < class TArr >
void 
SplitEdgeCallback< TArr >::
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

template < class TArr >
void
SplitEdgeCallback< TArr >::
reset( )
{
    this->hasFirstPoint = false;
    this->segmentGuide.setLine( 0, 0, 0, 0 );
    emit modelChanged( );
}

template < class TArr >
void
SplitEdgeCallback< TArr >::
slotModelChanged( )
{
    this->segmentGuide.update( );
}

template < class TArr >
void 
SplitEdgeCallback< TArr >::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{
    Point_2 clickedPoint = this->snapPoint( event );
    if ( ! this->hasFirstPoint )
    {
        this->p1 = clickedPoint;
        this->hasFirstPoint = true;
    }
    else
    {
        this->p2 = clickedPoint;
        X_monotone_curve_2 splitCurve =
            this->construct_x_monotone_curve_2( this->p1, this->p2 );
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
                if ( this->areEqual( hei->source( )->point( ), splitPoint ) ||
                    this->areEqual( hei->target( )->point( ), splitPoint ) )
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

template < class TArr >
void 
SplitEdgeCallback< TArr >::
mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{
    Point_2 clickedPoint = this->snapPoint( event );
    if ( this->hasFirstPoint )
    { // provide visual feedback for where the split line is
        Point_2 currentPoint = clickedPoint;
        Segment_2 currentSegment( this->p1, currentPoint );
        QLineF qSegment = this->convert( currentSegment );
        this->segmentGuide.setLine( qSegment );
        emit modelChanged( );
    }
}

template < class TArr >
typename SplitEdgeCallback< TArr >::Point_2
SplitEdgeCallback< TArr >::
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
