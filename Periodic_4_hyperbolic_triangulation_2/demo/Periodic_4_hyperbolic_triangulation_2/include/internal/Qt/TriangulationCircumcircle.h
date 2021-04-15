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

#ifndef CGAL_QT_TRIANGULATION_CIRCUMCIRCLE_H
#define CGAL_QT_TRIANGULATION_CIRCUMCIRCLE_H

#include <QGraphicsSceneMouseEvent>
#include <QGraphicsScene>
#include <QGraphicsEllipseItem>
#include <QEvent>
#include <QPen>

#include <CGAL/Qt/GraphicsViewInput.h>

namespace CGAL {
namespace Qt {

template <typename DT>
class TriangulationCircumcircle : public GraphicsViewInput
{
public:
  TriangulationCircumcircle(QGraphicsScene* s, DT  * dt_, QObject* parent);
  ~TriangulationCircumcircle();

  void setPen(const QPen& pen);

  void show();
  void hide();

protected:

  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

private:

  DT * dt;
  typename DT::Vertex_handle hint;
  typename DT::Face_handle fh;
  QGraphicsScene *scene_;
  QGraphicsEllipseItem* circle;
};


template <typename T>
TriangulationCircumcircle<T>::TriangulationCircumcircle(QGraphicsScene* s,
                                                              T * dt_,
                                                              QObject* parent)
  :  GraphicsViewInput(parent), dt(dt_), scene_(s)
{
  hint = typename T::Vertex_handle();
  circle = new QGraphicsEllipseItem();
  circle->hide();
  scene_->addItem(circle);
}


template <typename T>
TriangulationCircumcircle<T>::~TriangulationCircumcircle()
{
}


template <typename T>
void
TriangulationCircumcircle<T>::setPen(const QPen& pen)
{
  circle->setPen(pen);
}


template <typename T>
void
TriangulationCircumcircle<T>::show()
{
  circle->show();
}


template <typename T>
void
TriangulationCircumcircle<T>::hide()
{
  circle->hide();
}

template <typename T>
void
TriangulationCircumcircle<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  if(dt->dimension() != 2){
    circle->hide();
    return;
  }
  typename T::Point p = typename T::Point(event->scenePos().x(), event->scenePos().y());

  fh = dt->hyperbolic_locate(p);
  if(fh != typename T::Face_handle()) {
    hint = fh->vertex(0);
    typename T::Point p0, p1, p2;
    p0 = dt->geom_traits().construct_hyperbolic_point_2_object()(fh->vertex(0)->point(), fh->translation(0));
    p1 = dt->geom_traits().construct_hyperbolic_point_2_object()(fh->vertex(1)->point(), fh->translation(1));
    p2 = dt->geom_traits().construct_hyperbolic_point_2_object()(fh->vertex(2)->point(), fh->translation(2));

    typename T::Geom_traits::Circle_2 c(p0, p1, p2);
    CGAL::Bbox_2 bb = c.bbox();
    circle->setRect(bb.xmin(), bb.ymin(), bb.xmax()-bb.xmin(), bb.ymax()-bb.ymin());
    circle->show();
  } else {
    circle->hide();
  }
}


template <typename T>
bool
TriangulationCircumcircle<T>::eventFilter(QObject *obj, QEvent *event)
{
  if(event->type() == QEvent::GraphicsSceneMouseMove) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseMoveEvent(mouseEvent);
    return false; // don't consume the event
  } else{
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_TRIANGULATION_CIRCUMCIRCLE_H
