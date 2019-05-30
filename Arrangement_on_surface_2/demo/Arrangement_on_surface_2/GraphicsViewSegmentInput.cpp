// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "GraphicsViewSegmentInput.h"

#include <QGraphicsView>
#include <QEvent>

namespace CGAL {
namespace Qt {

GraphicsViewSegmentInputBase::
GraphicsViewSegmentInputBase( QObject* parent ):
  GraphicsViewInput( parent ),
  scene( NULL ),
  snappingEnabled( false ),
  snapToGridEnabled( false )
{ }

//! Setter function - sets the current scene
/*!
  \param scene_ A QGraphicsScene pointer to the class
*/
void GraphicsViewSegmentInputBase::setScene( QGraphicsScene* scene_ )
{
  this->scene = scene_;
}

//! Getter function - returns the current graph scene
/*!
  \return a QGraphicsScene object
*/
QGraphicsScene* GraphicsViewSegmentInputBase::getScene( ) const
{
  return this->scene;
}

//! It specifies the movement of the mouse on the graph scene
/*!
    \param QGraphicsSceneMouseEvent
*/
void GraphicsViewSegmentInputBase::
mouseMoveEvent(QGraphicsSceneMouseEvent* /* event */)
{ }

//! It specifies the mouse actions like drag,click
/*!
    \param QGraphicsSceneMouseEvent
*/
void GraphicsViewSegmentInputBase::
mousePressEvent(QGraphicsSceneMouseEvent* /* event */)
{
  // std::cout << "GraphicsViewSegmentInputBase::mousePressEvent" << std::endl;
}

//! determining the action of a mouse
/*!
  \param obj a QObject pointer to the class
  \param event a QEvent pointer to the class
  \return a boolean value 
*/
bool GraphicsViewSegmentInputBase::eventFilter( QObject* obj, QEvent* event )
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

//! Setting the snapping option 
/*!
  \param b a boolean value to toggle the current state
*/
void GraphicsViewSegmentInputBase::setSnappingEnabled( bool b )
{
  this->snappingEnabled = b;
}

//! Setting the snapping option to enabling grid
/*!
  \param b a boolean value to toggle the current state
*/
void GraphicsViewSegmentInputBase::setSnapToGridEnabled( bool b )
{
  this->snapToGridEnabled = b;
}

//! sets a new viewport box(rectangle) if not previously made
/*!
  \return A QRectF object
*/
QRectF GraphicsViewSegmentInputBase::viewportRect( ) const
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

} // namespace Qt
} // namespace CGAL
