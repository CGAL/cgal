#ifndef CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QEvent>
#include <QGraphicsLineItem>
#include <QGraphicsSceneMouseEvent>
#include <iostream>
#include "Callback.h"

namespace CGAL {
namespace Qt {

class GraphicsViewSegmentInputBase : public GraphicsViewInput
{
public:
    virtual void setScene( QGraphicsScene* scene_ );
    virtual QGraphicsScene* getScene( ) const;

protected:
    GraphicsViewSegmentInputBase( QObject* parent_ );
    virtual void mouseMoveEvent( QGraphicsSceneMouseEvent* event );
    virtual void mousePressEvent( QGraphicsSceneMouseEvent* event );
    virtual bool eventFilter( QObject* obj, QEvent* event );

    QGraphicsScene* scene;
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
void
GraphicsViewSegmentInput< K >::
mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{
    if ( this->second )
    {
        QPointF pt2 = event->scenePos( );
        QPointF pt = this->convert( this->p1 );
        this->segmentGuide.setLine( pt.x( ), pt.y( ), pt2.x( ), pt2.y( ) );
    }
}

template < class K >
void
GraphicsViewSegmentInput< K >::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{
    if ( !this->second )
    {
        this->second = 1;
        QPointF pt = event->scenePos( );
        this->p1 = this->convert( pt );
        this->segmentGuide.setLine( pt.x( ), pt.y( ), pt.x( ), pt.y( ) );
        if ( this->scene != NULL )
        {
            this->scene->addItem( &( this->segmentGuide ) );
        }
    }
    else
    {
        this->second = 0;
        this->p2 = this->convert( event->scenePos( ) );
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
