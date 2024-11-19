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

#ifndef CGAL_QT_GRAPHICS_VIEW_POLYLINE_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_POLYLINE_INPUT_H

#include <CGAL/license/GraphicsView.h>


#include <CGAL/auto_link/Qt.h>
#include <CGAL/export/Qt.h>

#include <QPolygonF>
#include <QPointF>
#include <QKeyEvent>

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsLineItem>

class QGraphicsScene;
class QGraphicsSceneMouseEvent;
class QGraphicsItem;
class QGraphicsPathItem;
class QGraphicsLineItem;
class QKeyEvent;
class QEvent;
class QObject;

namespace CGAL {
namespace Qt {

class CGAL_QT_EXPORT GraphicsViewPolylineInput_non_templated_base : public GraphicsViewInput
{
public:
  void setNumberOfVertices(int n)
  {
    n_ = n;
  }

  bool eventFilter(QObject *obj, QEvent *event);

  void setZValue(int v)
    {
      z = v;
    }

protected:
  // protected constructor
  GraphicsViewPolylineInput_non_templated_base(QObject* parent,
                                     QGraphicsScene* s,
                                     int n = 0,
                                     bool closed = true);


  // mousePressEvent returns true iff the event is consumed
  bool mousePressEvent(QGraphicsSceneMouseEvent *event);

  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);

  // keyPressEvent returns true iff the event is consumed
  bool keyPressEvent(QKeyEvent *event);

  void rubberbands(const QPointF& p);

  virtual void generate_polygon() = 0;

protected:
  QPolygonF polygon;
  bool closed_;

private:
  QGraphicsPathItem *path_item;
  QGraphicsLineItem *b, *e;
  int n_;
  QPointF sp;
  QGraphicsScene *scene_;
  int z;
}; // end class GraphicsViewPolylineInput_non_templated_base

template <typename K>
class GraphicsViewPolylineInput : public GraphicsViewPolylineInput_non_templated_base
{
public:
  GraphicsViewPolylineInput(QObject* parent, QGraphicsScene* s, int n = 0, bool closed = true)
    : GraphicsViewPolylineInput_non_templated_base(parent, s, n, closed)
  {
  }

protected:
  void generate_polygon() {
    std::list<typename K::Point_2> points;
    Converter<K> convert;
    convert(points, this->polygon);
    if(closed_ && points.size()>2){
      points.push_back(points.front());
    }
    Q_EMIT( generate(CGAL::make_object(points)));
  }
}; // end class GraphicsViewPolylineInput

} // namespace Qt
} // namespace CGAL

#ifdef CGAL_HEADER_ONLY
#include <CGAL/Qt/GraphicsViewPolylineInput_impl.h>
#endif // CGAL_HEADER_ONLY

#endif // CGAL_QT_GRAPHICS_VIEW_POLYLINE_INPUT_H
