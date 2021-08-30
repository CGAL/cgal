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
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_POINT_INPUT
#define CGAL_QT_POINT_INPUT

#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsSceneMouseEvent>
#include <QEvent>
#include <list>



namespace CGAL {
namespace Qt {

template <typename DT>
class TriangulationPointInput : public GraphicsViewInput
{
public:
  typedef typename DT::Geom_traits K;
  typedef typename DT::Face_handle Face_handle;
  typedef typename DT::Point Point;

  TriangulationPointInput(QGraphicsScene* s, DT  * dt_, QObject* parent);

protected:
  void localize_and_insert_point(QPointF qt_point);

  void mousePressEvent(QGraphicsSceneMouseEvent *event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
  bool eventFilter(QObject *obj, QEvent *event);

  std::list<Face_handle> faces;
  std::list<QGraphicsPolygonItem*> qfaces;
  DT * dt;
  Converter<K> convert;
  QGraphicsScene *scene_;
  Point p;
};


template <typename T>
TriangulationPointInput<T>::TriangulationPointInput(QGraphicsScene* s,
              T * dt_,
              QObject* parent)
  :  GraphicsViewInput(parent), dt(dt_), scene_(s)
{}




// TODO: Do this!

template <typename T>
void
TriangulationPointInput<T>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
  p = convert(event->scenePos());
}



template <typename T>
void
TriangulationPointInput<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent * /*event*/)
{
  faces.clear();
  for(std::list<QGraphicsPolygonItem*>::iterator it = qfaces.begin();
      it != qfaces.end();
      ++it){
    scene_->removeItem(*it);
    delete *it;
  }
  qfaces.clear();
  emit (generate(CGAL::make_object(p)));
}



template <typename T>
bool
TriangulationPointInput<T>::eventFilter(QObject *obj, QEvent *event)
{
  if(event->type() == QEvent::GraphicsSceneMousePress) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mousePressEvent(mouseEvent);
    return true;
  } else if(event->type() == QEvent::GraphicsSceneMouseRelease) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseReleaseEvent(mouseEvent);
    return true;
  } else{
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_POINT_INPUT
