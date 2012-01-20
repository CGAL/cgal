// Copyright (c) 2008  GeometryFactory Sarl (France).
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
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_GRAPHICS_VIEW_ISO_RECTANGLE_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_ISO_RECTANGLE_INPUT_H

#include <QGraphicsView>
#include <QRectF>
#include <QPointF>
#include <QGraphicsItem>
#include <QGraphicsRectItem> 
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>

#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/GraphicsViewInput.h>

#include <CGAL/array.h>

namespace CGAL {
namespace Qt {

template <typename K>
class GraphicsViewIsoRectangleInput : public GraphicsViewInput
{
public:
  GraphicsViewIsoRectangleInput(QObject *parent, QGraphicsScene* s); 
  ~GraphicsViewIsoRectangleInput();

protected:
    
  virtual void mousePressEvent(QGraphicsSceneMouseEvent *event);
  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  virtual void keyPressEvent(QKeyEvent *event);
  
  bool eventFilter(QObject *obj, QEvent *event);
  

  

private:

  typedef typename K::Point_2 Point_2;
  QPointF qp, qq, qr;
  Point_2 p, q, r;
  QGraphicsRectItem *rectItem;
  QPointF rect_first_point;
  QGraphicsScene *scene_;  
  Converter<K> convert;
};


template <typename K>
GraphicsViewIsoRectangleInput<K>::GraphicsViewIsoRectangleInput(QObject *parent, QGraphicsScene* s)
  : GraphicsViewInput(parent), rectItem(new QGraphicsRectItem), scene_(s)
{
  rectItem->setBrush(QBrush());
  scene_->addItem(rectItem);
  rectItem->hide();
  rectItem->setZValue(10000);
}

template <typename K>
GraphicsViewIsoRectangleInput<K>::~GraphicsViewIsoRectangleInput()
{
  delete rectItem;
}


template <typename K>
void 
GraphicsViewIsoRectangleInput<K>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{  
  if(event->modifiers()  & ::Qt::ShiftModifier){
    return;
  }
  if(event->button() == ::Qt::LeftButton) {
    if(rectItem->isVisible()) {
      // we have clicked twice
      emit generate(CGAL::make_object(convert(rectItem->rect())));
      rectItem->hide();
    } else { 
      // we enter a first point
      rect_first_point = event->scenePos();
      rectItem->setRect(QRectF(rect_first_point, rect_first_point));
      rectItem->show();
    }
  }
}


template <typename K>
void 
GraphicsViewIsoRectangleInput<K>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  // todo: only do this if no modifiers are pressed at the same time
  if(rectItem->isVisible()) {
    rectItem->setRect(QRectF(rect_first_point,
			     event->scenePos()));
  }
}


template <typename K>
void 
GraphicsViewIsoRectangleInput<K>::keyPressEvent ( QKeyEvent * event ) 
{
}



template <typename K>
bool 
GraphicsViewIsoRectangleInput<K>::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::GraphicsSceneMousePress) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mousePressEvent(mouseEvent);
    return true;
  } else if (event->type() == QEvent::GraphicsSceneMouseMove) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseMoveEvent(mouseEvent);
    return true;
  } else if (event->type() == QEvent::KeyPress) {
    QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
    keyPressEvent(keyEvent);
    return true;
  } else{
    // standard event processing
    return QObject::eventFilter(obj, event);
  }
} 

} // namespace Qt

} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_ISO_RECTANGLE_INPUT_H
