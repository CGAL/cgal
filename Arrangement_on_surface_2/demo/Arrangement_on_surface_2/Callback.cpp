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

#include "Callback.h"

#include <QEvent>
#include <QKeyEvent>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

namespace CGAL {
namespace Qt {

Callback::Callback( QObject* parent ) : QObject( parent ) { }

void Callback::reset( ) { }

bool Callback::eventFilter( QObject* object, QEvent* event )
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
    QKeyEvent* keyEvent = static_cast< QKeyEvent* >( event );
    this->keyPressEvent( keyEvent );
  }
  return QObject::eventFilter( object, event );
}

  void Callback::mousePressEvent(QGraphicsSceneMouseEvent* /* event */) { }

  void Callback::mouseMoveEvent(QGraphicsSceneMouseEvent* /* event */) { }

  void Callback::mouseReleaseEvent(QGraphicsSceneMouseEvent* /* event */) { }

  void Callback::keyPressEvent(QKeyEvent* /* event */)  { }

void Callback::slotModelChanged( ) { }

} // namespace Qt
} // namespace CGAL
