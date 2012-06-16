#include "Callback.h"
#include <QEvent>
#include <QKeyEvent>
#include <QGraphicsSceneMouseEvent>

namespace CGAL {
namespace Qt {

Callback::
Callback( QObject* parent ) 
    : QObject( parent )
{}

bool 
Callback::
eventFilter( QObject* object, QEvent* event )
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
        this->mouseReleaseEvent( mouseEvent );
    }
    else if ( event->type( ) == QEvent::KeyPress )
    {
        QKeyEvent* keyEvent =
            static_cast< QKeyEvent* >( event );
        this->keyPressEvent( keyEvent );
    }
    return QObject::eventFilter( object, event );
}

void
Callback::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{

}

void
Callback::
mouseMoveEvent( QGraphicsSceneMouseEvent* event )
{

}

void
Callback::
mouseReleaseEvent( QGraphicsSceneMouseEvent* event )
{

}

void
Callback::
keyPressEvent( QKeyEvent* event )
{

}

} // namespace Qt
} // namespace CGAL
