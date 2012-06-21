#include "Callback.h"
#include <QEvent>
#include <QKeyEvent>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

namespace CGAL {
namespace Qt {

Callback::
Callback( QObject* parent ):
    QObject( parent ),
    scene( NULL )
{}

void 
Callback::
setScene( QGraphicsScene* scene_ )
{
    this->scene = scene_;
}

void 
Callback::
reset( )
{

}

QGraphicsScene* 
Callback::
getScene( ) const
{
    return this->scene;
}

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

QRectF
Callback::
viewportRect( ) const
{
    QRectF res;
    if ( this->scene == NULL )
    {
        return res;
    }

    QList< QGraphicsView* > views = this->scene->views( );
    if ( views.size( ) == 0 )
    {
        return res;
    }
    // assumes the first view is the right one
    QGraphicsView* viewport = views.first( );
    QPointF p1 = viewport->mapToScene( 0, 0 );
    QPointF p2 = viewport->mapToScene( viewport->width( ), viewport->height( ) );
    res = QRectF( p1, p2 );

    return res;
}


void
Callback::
slotModelChanged( )
{

}

} // namespace Qt
} // namespace CGAL
