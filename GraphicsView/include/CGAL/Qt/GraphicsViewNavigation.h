// Copyright (c) 2008  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_GRAPHICS_VIEW_NAVIGATION_H
#define CGAL_QT_GRAPHICS_VIEW_NAVIGATION_H

#include <CGAL/license/GraphicsView.h>


#include <CGAL/auto_link/Qt.h>
#include <CGAL/export/Qt.h>

#include <QObject>
#include <QPointF>
#include <QString>
#include <QCursor>
#include <QRect>
#include <QRectF>

class QGraphicsView;
class QGraphicsScene;
class QEvent;
class QGraphicsRectItem;

namespace CGAL {
namespace Qt {

class CGAL_QT_EXPORT GraphicsViewNavigation: public QObject {

  Q_OBJECT

  Q_SIGNALS:
  void mouseCoordinates(QString);

public:
  GraphicsViewNavigation();
  ~GraphicsViewNavigation();

  bool eventFilter(QObject *obj, QEvent *event);

private:

  void scaleView(QGraphicsView*, qreal scaleFactor);
  void translateView(QGraphicsView*, int dx,  int dy);
  void drag_to(QGraphicsView*, QPoint new_pos);
  void display_parameters(QGraphicsView*);

  QGraphicsRectItem* rectItem;
  QPointF rect_first_point;
  bool dragging;
  QPointF dragging_start;
  QCursor cursor_backup;
};

} // namespace Qt
} // namespace CGAL

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Qt/GraphicsViewNavigation_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_QT_GRAPHICS_VIEW_NAVIGATION_H
