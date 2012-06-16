#ifndef DELETE_CURVE_CALLBACK_H
#define DELETE_CURVE_CALLBACK_H
#include "Callback.h"
#include <QEvent>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/CurveGraphicsItem.h>
#include <CGAL/Arrangement_with_history_2.h>
#include "Utils.h"

/**
Handles deletion of arrangement curves selected from the scene.

The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < class TArr >
class DeleteCurveCallback : public CGAL::Qt::Callback
{
public:
    typedef typename TArr::Halfedge_handle Halfedge_handle;
    typedef typename TArr::Halfedge_iterator Halfedge_iterator;
    typedef typename TArr::Geometry_traits_2 Traits;
    typedef typename TArr::Curve_handle Curve_handle;
    typedef typename TArr::Originating_curve_iterator Originating_curve_iterator;
    typedef typename TArr::Induced_edge_iterator Induced_edge_iterator;
    typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
    typedef typename Traits::Kernel Kernel;
    typedef typename Kernel::Point_2 Point;
    typedef typename Kernel::Segment_2 Segment;

    DeleteCurveCallback( TArr* arr_, QObject* parent_ );
    void setScene( QGraphicsScene* scene_ );
    QGraphicsScene* getScene( ) const;
    void reset( );

protected:
    void mousePressEvent( QGraphicsSceneMouseEvent *event );
    void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
    void highlightNearestCurve( QGraphicsSceneMouseEvent *event );

    Compute_squared_distance_2< Traits > squaredDistance;
    CGAL::Qt::Converter< Kernel > convert;
    QGraphicsScene* scene;
    CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurve;
    TArr* arr;
    Halfedge_handle removableHalfedge;
}; // class DeleteCurveCallback


template < class TArr >
DeleteCurveCallback< TArr >::
DeleteCurveCallback( TArr* arr_, QObject* parent_ ):
    CGAL::Qt::Callback( parent_ ),
    arr( arr_ ),
    scene( NULL ),
    highlightedCurve( new CGAL::Qt::CurveGraphicsItem< Traits >( ) )
{
    QObject::connect( this, SIGNAL( modelChanged( ) ),
        this->highlightedCurve, SLOT( modelChanged( ) ) );
}

template < class TArr >
void 
DeleteCurveCallback< TArr >::
setScene( QGraphicsScene* scene_ )
{
    this->scene = scene_;
    if ( this->scene )
    {
        this->scene->addItem( this->highlightedCurve );
    }
}

template < class TArr >
QGraphicsScene* 
DeleteCurveCallback< TArr >::
getScene( ) const
{
    return this->scene;
}

template < class TArr >
void
DeleteCurveCallback< TArr >::
reset( )
{
    this->highlightedCurve->clear( );
    this->removableHalfedge = Halfedge_handle( );
    emit modelChanged( );
}

template < class TArr >
void 
DeleteCurveCallback< TArr >::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{
    if ( this->removableHalfedge == Halfedge_handle( ) )
    {
        return;
    }

    Originating_curve_iterator it = this->arr->originating_curves_begin( this->removableHalfedge );
    Originating_curve_iterator it_end = this->arr->originating_curves_end( this->removableHalfedge );
    while ( it != it_end )
    {
        Originating_curve_iterator temp = it;
        ++temp;
        CGAL::remove_curve( *(this->arr), it );
        it = temp;
    }

    this->reset( );
}

template < class TArr >
void 
DeleteCurveCallback< TArr >::
mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{
    this->highlightNearestCurve( event );
}

template < class TArr >
void 
DeleteCurveCallback< TArr >::
highlightNearestCurve( QGraphicsSceneMouseEvent* event )
{
    // find the nearest curve to the cursor to be the new highlighted curve
    Point p = this->convert( event->scenePos( ) );
    bool isFirst = true;
    double minDist = 0.0;
    Halfedge_iterator nearestHei;

    for ( Halfedge_iterator hei = this->arr->halfedges_begin( );
        hei != this->arr->halfedges_end( );
        ++hei )
    {
        X_monotone_curve_2 curve = hei->curve( );
        double dist = CGAL::to_double( this->squaredDistance( p, curve ) );
        std::cout << dist << std::endl;
        if ( isFirst || dist < minDist )
        {
            isFirst = false;
            minDist = dist;
            nearestHei = hei;
        }
    }

    // now 'removableHalfedge' holds the closest halfedge to the point of the mouse
    this->removableHalfedge = nearestHei;
    if ( isFirst )
    {
        return;
    }

    // create a curve graphics item and add it to the scene
    this->highlightedCurve->clear( );
    Originating_curve_iterator ocit, temp;
    ocit = this->arr->originating_curves_begin( nearestHei );
    while ( ocit != this->arr->originating_curves_end( nearestHei ) )
    {
        temp = ocit;
        ++temp;

        Curve_handle ch = ocit;
        Induced_edge_iterator itr;
        for ( itr = this->arr->induced_edges_begin( ch );
                itr != this->arr->induced_edges_end( ch );
                ++itr )
        {
            X_monotone_curve_2 curve = (*itr)->curve( );
            this->highlightedCurve->insert( curve );
        }
        ocit = temp;
    }

    emit modelChanged( );
}

#endif // DELETE_CURVE_CALLBACK_H
