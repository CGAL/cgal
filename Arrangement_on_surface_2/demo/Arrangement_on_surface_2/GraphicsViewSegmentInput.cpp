#include "GraphicsViewSegmentInput.h"

namespace CGAL {
namespace Qt {

GraphicsViewSegmentInputBase::
GraphicsViewSegmentInputBase( QObject* parent ):
    GraphicsViewInput( parent ),
    scene( NULL )
{ }

void 
GraphicsViewSegmentInputBase::
setScene( QGraphicsScene* scene_ )
{
    this->scene = scene_;
}

QGraphicsScene* 
GraphicsViewSegmentInputBase::
getScene( ) const
{
    return this->scene;
}

void 
GraphicsViewSegmentInputBase::
mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{ }

void 
GraphicsViewSegmentInputBase::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{ }

bool 
GraphicsViewSegmentInputBase::
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

    return QObject::eventFilter( obj, event );
}

} // namespace Qt
} // namespace CGAL
