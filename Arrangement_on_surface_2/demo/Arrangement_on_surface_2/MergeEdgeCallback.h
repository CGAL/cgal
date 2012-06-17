#ifndef MERGE_EDGE_CALLBACK_H
#define MERGE_EDGE_CALLBACK_H
#include "Callback.h"
#include <QEvent>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/CurveGraphicsItem.h>
#include <CGAL/Arrangement_with_history_2.h>
#include "Utils.h"

/**
Handles merging of arrangement curves selected from the scene.

The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < class TArr >
class MergeEdgeCallback : public CGAL::Qt::Callback
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
    typedef typename Traits::Kernel Kernel;
    typedef typename Kernel::Point_2 Point;
    typedef typename Kernel::Segment_2 Segment;

    MergeEdgeCallback( TArr* arr_, QObject* parent_ );
    void setScene( QGraphicsScene* scene_ );
    QGraphicsScene* getScene( ) const;
    void reset( );

protected:
    void mousePressEvent( QGraphicsSceneMouseEvent *event );
    void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
    Halfedge_handle getNearestMergeableCurve( QGraphicsSceneMouseEvent *event );
    Halfedge_handle getNearestMergeableCurve( Halfedge_handle h, QGraphicsSceneMouseEvent *event );

    Compute_squared_distance_2< Traits > squaredDistance;
    CGAL::Qt::Converter< Kernel > convert;
    QGraphicsScene* scene;
    CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurve;
    CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurve2;
    TArr* arr;
    Halfedge_handle mergeableHalfedge;
    bool isFirst;
}; // class MergeEdgeCallback


template < class TArr >
MergeEdgeCallback< TArr >::
MergeEdgeCallback( TArr* arr_, QObject* parent_ ):
    CGAL::Qt::Callback( parent_ ),
    arr( arr_ ),
    scene( NULL ),
    highlightedCurve( new CGAL::Qt::CurveGraphicsItem< Traits >( ) ),
    highlightedCurve2( new CGAL::Qt::CurveGraphicsItem< Traits >( ) ),
    isFirst( true )
{
    QObject::connect( this, SIGNAL( modelChanged( ) ),
        this->highlightedCurve, SLOT( modelChanged( ) ) );
    QObject::connect( this, SIGNAL( modelChanged( ) ),
        this->highlightedCurve2, SLOT( modelChanged( ) ) );
}

template < class TArr >
void 
MergeEdgeCallback< TArr >::
setScene( QGraphicsScene* scene_ )
{
    this->scene = scene_;
    if ( this->scene )
    {
        this->scene->addItem( this->highlightedCurve );
        this->scene->addItem( this->highlightedCurve2 );
    }
}

template < class TArr >
QGraphicsScene* 
MergeEdgeCallback< TArr >::
getScene( ) const
{
    return this->scene;
}

template < class TArr >
void
MergeEdgeCallback< TArr >::
reset( )
{
    this->isFirst = true;
    this->highlightedCurve->clear( );
    this->highlightedCurve2->clear( );
    this->mergeableHalfedge = Halfedge_handle( );
    emit modelChanged( );
}

template < class TArr >
void 
MergeEdgeCallback< TArr >::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{
    if ( this->isFirst )
    { // save the first edge if mergeable
        Halfedge_handle halfedge = this->getNearestMergeableCurve( event );
        if ( halfedge == Halfedge_handle( ) )
        {
            return;
        }
        this->isFirst = false;
        this->mergeableHalfedge = halfedge;
    }
    else
    {
        Halfedge_handle nextHalfedge = 
            this->getNearestMergeableCurve( this->mergeableHalfedge, event );
        this->arr->merge_edge( this->mergeableHalfedge, nextHalfedge );
        this->reset( );
    }

    emit modelChanged( );
}

template < class TArr >
void 
MergeEdgeCallback< TArr >::
mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{
    if ( this->isFirst )
    {
        Halfedge_handle halfedge = this->getNearestMergeableCurve( event );
        if ( halfedge == Halfedge_handle( ) )
        {
            return;
        }
        this->highlightedCurve->clear( );
        this->highlightedCurve->insert( halfedge->curve( ) );
        emit modelChanged( );
    }
    else
    {
        Halfedge_handle nextHalfedge = 
            this->getNearestMergeableCurve( this->mergeableHalfedge, event );
        this->highlightedCurve2->clear( );
        this->highlightedCurve2->insert( nextHalfedge->curve( ) );
        emit modelChanged( );
    }
}

template < class TArr >
typename MergeEdgeCallback< TArr >::Halfedge_handle
MergeEdgeCallback< TArr >::
getNearestMergeableCurve( QGraphicsSceneMouseEvent* event )
{
    // find the nearest curve to the cursor that is adjacent to a curve that
    // can be merged with it
    Point p = this->convert( event->scenePos( ) );
    double minDist = 0.0;
    bool noneFound = true;
    Halfedge_iterator nearestHei;

    for ( Halfedge_iterator hei = this->arr->halfedges_begin( );
        hei != this->arr->halfedges_end( );
        ++hei )
    {
        Vertex_iterator source = hei->source( );
        Vertex_iterator target = hei->target( );
        if ( source->degree( ) != 2 && target->degree( ) != 2 )
        { // then this halfedge has no mergeable neighbors
            continue;
        }

        X_monotone_curve_2 curve = hei->curve( );
        double dist = CGAL::to_double( this->squaredDistance( p, curve ) );
        if ( noneFound || dist < minDist )
        {
            noneFound = false;
            minDist = dist;
            nearestHei = hei;
        }
    }

    if ( noneFound )
    { // then we did not find a mergeable halfedge
        return Halfedge_handle( );
    }
    return nearestHei;
}

template < class TArr >
typename MergeEdgeCallback< TArr >::Halfedge_handle
MergeEdgeCallback< TArr >::
getNearestMergeableCurve( Halfedge_handle h, QGraphicsSceneMouseEvent* event )
{
    // find the nearest curve to the cursor that is adjacent to a curve that
    // can be merged with it
    Point p = this->convert( event->scenePos( ) );
    Halfedge_handle h1 = h->prev( );
    Halfedge_handle h2 = h->next( );
    Vertex_iterator source = h->source( );
    Vertex_iterator target = h->target( );
    if ( source->degree( ) != 2 && target->degree( ) != 2 )
    {
        return Halfedge_handle( );
    }
    else if ( source->degree( ) != 2 )
    {
        return h2;
    }
    else if ( target->degree( ) != 2 )
    {
        return h1;
    }
    else
    {
        X_monotone_curve_2 c1 = h1->curve( );
        X_monotone_curve_2 c2 = h2->curve( );
        double d1 = CGAL::to_double( this->squaredDistance( p, c1 ) );
        double d2 = CGAL::to_double( this->squaredDistance( p, c2 ) );
        if ( d1 < d2 )
        {
            return h1;
        }
        else
        {
            return h2;
        }
    }
}
#endif // MERGE_EDGE_CALLBACK_H
