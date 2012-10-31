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

#ifndef CGAL_QT_CALLBACK_H
#define CGAL_QT_CALLBACK_H

#include <QObject>

#include "Utils.h"

class QRectF;
class QEvent;
class QKeyEvent;
class QGraphicsScene;
class QGraphicsSceneMouseEvent;

namespace CGAL {
namespace Qt {

class Callback : public QObject, public QGraphicsSceneMixin
{
Q_OBJECT

public:
  Callback( QObject* parent );
  virtual void reset( );

public slots:
  virtual void slotModelChanged( );

signals:
  void modelChanged( );

protected:
  virtual bool eventFilter( QObject* object, QEvent* event );
  virtual void mousePressEvent( QGraphicsSceneMouseEvent* event );
  virtual void mouseMoveEvent( QGraphicsSceneMouseEvent* event );
  virtual void mouseReleaseEvent( QGraphicsSceneMouseEvent* event );
  virtual void keyPressEvent( QKeyEvent* event );
};

} // namespace Qt
} // namespace CGAL
#endif // CGAL_QT_CALLBACK_H
