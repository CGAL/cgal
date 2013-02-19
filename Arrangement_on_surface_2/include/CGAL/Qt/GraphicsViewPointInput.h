// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$

#ifndef CGAL_QT_GRAPHICS_VIEW_POINT_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_POINT_INPUT_H
#include <CGAL/Qt/GraphicsViewInput.h>
#include <QEvent>
#include <QGraphicsSceneMouseEvent>
#include <iostream>
namespace CGAL {
namespace Qt {
template < class K >
class GraphicsViewPointInput: public GraphicsViewInput
{
public:
    typedef typename K::Point_2 Point;

    GraphicsViewPointInput( QObject* parent );

protected:
    void mousePressEvent( QGraphicsSceneMouseEvent* event );
    void mouseReleaseEvent( QGraphicsSceneMouseEvent* event );
    bool eventFilter( QObject* obj, QEvent* event );

    Converter< K > convert;
    Point point;

}; // class GraphicsViewPointInput

template < class K >
GraphicsViewPointInput< K >::
GraphicsViewPointInput( QObject* parent ):
    GraphicsViewInput( parent )
{

}

template < class K >
void
GraphicsViewPointInput< K >::
mousePressEvent( QGraphicsSceneMouseEvent* event )
{
    std::cout << event->pos( ).x( ) << std::endl;
    this->point = this->convert( event->scenePos( ) );
}

template < class K >
void
GraphicsViewPointInput< K >::
mouseReleaseEvent( QGraphicsSceneMouseEvent* event )
{
    emit generate( CGAL::make_object( this->point ) );
}

template < class K >
bool
GraphicsViewPointInput< K >::
eventFilter( QObject* obj, QEvent* event )
{
    if ( event->type( ) == QEvent::GraphicsSceneMousePress )
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
    return QObject::eventFilter( obj, event );
}

} // namespace Qt
} // namespace CGAL
#endif // CGAL_QT_GRAPHICS_VIEW_POINT_INPUT_H
