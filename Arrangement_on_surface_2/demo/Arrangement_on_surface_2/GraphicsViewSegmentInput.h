#ifndef CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
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

class GraphicsViewSegmentInputBase : public GraphicsViewInput, public ISnappable
{
public:
    virtual void setScene( QGraphicsScene* scene_ );
    QGraphicsScene* getScene( ) const;

    void setSnappingEnabled( bool b );
    void setSnapToGridEnabled( bool b );

protected:
    GraphicsViewSegmentInputBase( QObject* parent_ );
    virtual void mouseMoveEvent( QGraphicsSceneMouseEvent* event );
    virtual void mousePressEvent( QGraphicsSceneMouseEvent* event );
    virtual bool eventFilter( QObject* obj, QEvent* event );

    QRectF viewportRect( ) const;

    QGraphicsScene* scene;
    bool snappingEnabled;
    bool snapToGridEnabled;

}; // class GraphicsViewSegmentInputBase

template < class K >
class GraphicsViewSegmentInput: public GraphicsViewSegmentInputBase
{
public:
    typedef typename K::Point_2 Point;
    typedef typename K::Segment_2 Segment;

    GraphicsViewSegmentInput( QObject* parent );

protected:
    void mouseMoveEvent( QGraphicsSceneMouseEvent* event );
    void mousePressEvent( QGraphicsSceneMouseEvent* event );

    // override this to snap to the points you like
    virtual Point snapPoint( QGraphicsSceneMouseEvent* event );

    Converter< K > convert;
    Point p1;
    Point p2;
    bool second;

    QGraphicsLineItem segmentGuide;
}; // class GraphicsViewSegmentInput


template < class K >
GraphicsViewSegmentInput< K >::
GraphicsViewSegmentInput( QObject* parent ):
    GraphicsViewSegmentInputBase( parent ),
    second( false )
{ }

template < class K >
typename K::Point_2
GraphicsViewSegmentInput< K >::
snapPoint( QGraphicsSceneMouseEvent* event )
{
    Point clickedPoint = this->convert( event->scenePos( ) );
    return clickedPoint;
}

template < class K >
void
GraphicsViewSegmentInput< K >::
mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{
    if ( this->second )
    {
        Point clickedPoint = this->snapPoint( event );
        Segment segment( this->p1, clickedPoint );
        QLineF qSegment = this->convert( segment );
        this->segmentGuide.setLine( qSegment );
    }
}

template < class K >
void
GraphicsViewSegmentInput< K >::
mousePressEvent( QGraphicsSceneMouseEvent* event )
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
        Segment res( this->p1, this->p2 );
        emit generate( CGAL::make_object( res ) );
    }
}

} // namespace Qt
} // namespace CGAL
#endif // CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
