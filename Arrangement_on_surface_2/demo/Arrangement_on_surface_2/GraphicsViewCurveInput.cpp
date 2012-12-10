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
// $URL: $
// $Id: $
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "GraphicsViewCurveInput.h"

#include <QGraphicsView>

namespace CGAL {
namespace Qt {

void GraphicsViewCurveInputBase::setScene( QGraphicsScene* scene_ )
{
  this->QGraphicsSceneMixin::setScene( scene_ );
  if ( this->scene != NULL )
  {
    this->scene->addItem( &this->pointsGraphicsItem );
  }
}

#if 0
QGraphicsScene* GraphicsViewCurveInputBase::getScene( ) const
{
  return this->scene;
}
#endif

void GraphicsViewCurveInputBase::setSnappingEnabled( bool b )
{
  this->snappingEnabled = b;
}

void GraphicsViewCurveInputBase::setSnapToGridEnabled( bool b )
{
  this->snapToGridEnabled = b;
}

/*! Constructor */
GraphicsViewCurveInputBase::GraphicsViewCurveInputBase( QObject* parent ) :
    GraphicsViewInput( parent ),
//    scene( NULL ),
    snappingEnabled( false ),
    snapToGridEnabled( false ),
    color( ::Qt::blue )
{
    this->pointsGraphicsItem.setZValue( 100 );
    this->pointsGraphicsItem.setColor( this->color );
}

void
GraphicsViewCurveInputBase::mouseMoveEvent(QGraphicsSceneMouseEvent* /* event */)
{ }

void GraphicsViewCurveInputBase::
mousePressEvent(QGraphicsSceneMouseEvent* /* event */)
{ 
  // std::cout << "GraphicsViewCurveInputBase::mousePressEvent" << std::endl;
}

bool GraphicsViewCurveInputBase::eventFilter( QObject* obj, QEvent* event )
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

void GraphicsViewCurveInputBase::setColor( QColor c )
{
  this->color = c;
  this->pointsGraphicsItem.setColor( this->color );
}

QColor GraphicsViewCurveInputBase::getColor( ) const
{
  return this->color;
}

#if 0
QRectF GraphicsViewCurveInputBase::viewportRect( ) const
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
#endif

} // namespace Qt
} // namespace CGAL
