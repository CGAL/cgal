// Copyright (c) 2008  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//

#ifndef CGAL_QT_TRIANGULATION_CONFLICT_ZONE
#define CGAL_QT_TRIANGULATION_CONFLICT_ZONE

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/TriangulationGraphicsItem.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <list>

namespace CGAL {
namespace Qt {

template <typename DT>
class TriangulationConflictZone : public GraphicsViewInput
{
public:
  typedef typename DT::Geom_traits  K;
  typedef typename DT::Face_handle  Face_handle;
  typedef typename DT::Point        Point;
  typedef typename DT::Hyperbolic_translation       Hyperbolic_translation;
  typedef typename DT::Hyperbolic_triangle     Triangle;

  TriangulationConflictZone(QGraphicsScene* s, DT  * dt_, QObject* parent);

protected:

  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

  Point p;

  DT * dt;
  Converter<K> convert;
  QGraphicsScene *scene_;
};


template <typename T>
TriangulationConflictZone<T>::TriangulationConflictZone(QGraphicsScene* s,
              T * dt_,
              QObject* parent)
  :  GraphicsViewInput(parent), dt(dt_), scene_(s)
{}



template <typename T>
void
TriangulationConflictZone<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  p = convert(event->scenePos());
  emit(generate(CGAL::make_object(p)));
}


template <typename T>
bool
TriangulationConflictZone<T>::eventFilter(QObject *obj, QEvent *event)
{
  if(event->type() == QEvent::GraphicsSceneMouseMove) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseMoveEvent(mouseEvent);
    return false; // don't consume the event
  } else {
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_TRIANGULATION_CONFLICT_ZONE
