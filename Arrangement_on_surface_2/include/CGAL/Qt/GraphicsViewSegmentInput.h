#ifndef CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
#include <CGAL/Qt/GraphicsViewInput.h>
#include <QEvent>
#include <QGraphicsSceneMouseEvent>
#include <iostream>
namespace CGAL {
namespace Qt {
template < class K >
class GraphicsViewSegmentInput: public GraphicsViewInput
{
public:
    typedef typename K::Point_2 Point;
    typedef typename K::Segment_2 Segment;

    GraphicsViewSegmentInput( QObject* parent );

    void setScene( QGraphicsScene* scene_ );

protected:
    void mouseMoveEvent( QGraphicsSceneMouseEvent* event );
    void mousePressEvent( QGraphicsSceneMouseEvent* event );
    void mouseReleaseEvent( QGraphicsSceneMouseEvent* event );
    bool eventFilter( QObject* obj, QEvent* event );

    Converter< K > convert;
    Point p1;
    Point p2;
    bool second;

    QGraphicsLineItem segmentGuide;
    QGraphicsScene* scene;
}; // class GraphicsViewSegmentInput

template < class K >
GraphicsViewSegmentInput< K >::
GraphicsViewSegmentInput( QObject* parent ):
    GraphicsViewInput( parent ),
    second( 0 ),
    scene( NULL )
{ }

template < class K >
void
GraphicsViewSegmentInput< K >::
setScene( QGraphicsScene* scene_ )
{
    this->scene = scene_;
}

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
        //this->p2 = this->convert( event->scenePos( ) );
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

template < class K >
void
GraphicsViewSegmentInput< K >::
mouseReleaseEvent( QGraphicsSceneMouseEvent* event )
{
    if ( !this->second )
    {

    }
    else
    {
    }
}

template < class K >
bool
GraphicsViewSegmentInput< K >::
eventFilter( QObject* obj, QEvent* event )
{
    if ( event->type( ) == QEvent::GraphicsSceneMouseMove )
    {
        QGraphicsSceneMouseEvent* mouseEvent =
            static_cast< QGraphicsSceneMouseEvent* >( event );
        this->mouseMoveEvent( mouseEvent );
    }
    else if ( event->type( ) == QEvent::GraphicsSceneMousePress )
    {
        QGraphicsSceneMouseEvent* mouseEvent =
            static_cast< QGraphicsSceneMouseEvent* >( event );
        this->mousePressEvent( mouseEvent );
    }
    else if ( event->type( ) == QEvent::GraphicsSceneMouseRelease )
    {
        QGraphicsSceneMouseEvent* mouseEvent =
            static_cast< QGraphicsSceneMouseEvent* >( event );
        //this->mouseReleaseEvent( mouseEvent );
    }
    return QObject::eventFilter( obj, event );
}

} // namespace Qt
} // namespace CGAL
#endif // CGAL_QT_GRAPHICS_VIEW_SEGMENT_INPUT_H
