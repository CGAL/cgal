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

#ifndef CGAL_QT_GRAPHICS_VIEW_CIRCULAR_ARC_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_CIRCULAR_ARC_INPUT_H

#include <QGraphicsView>
#include <QRectF>
#include <QPointF>
#include <QGraphicsItem>
#include <QGraphicsEllipseItem> 
#include <QGraphicsLineItem> 
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>

#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/CircularArcGraphicsItem.h>

#include <CGAL/array.h>

namespace CGAL {
namespace Qt {

template <typename K>
class GraphicsViewCircularArcInput : public GraphicsViewInput
{
public:
  GraphicsViewCircularArcInput(QObject *parent, QGraphicsScene* s); 
  ~GraphicsViewCircularArcInput();

protected:
    
  virtual void mousePressEvent(QGraphicsSceneMouseEvent *event);
  virtual void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
  virtual void keyPressEvent(QKeyEvent *event);
  
  bool eventFilter(QObject *obj, QEvent *event);
  

  

private:
  typedef typename K::Circular_arc_2 Circular_arc_2;
  typedef typename K::Line_arc_2 Line_arc_2;
  typedef typename K::Point_2 Point_2;
  typedef typename K::Segment_2 Segment_2;

  int count;
  QGraphicsLineItem *qline;
  CircularArcGraphicsItem<K> *qcarc;
  QPointF qp, qq, qr;
  Point_2 p, q, r, ap, aq, ar;
  QGraphicsScene *scene_;  
  Converter<K> convert;
};


template <typename K>
GraphicsViewCircularArcInput<K>::GraphicsViewCircularArcInput(QObject *parent, QGraphicsScene* s)
  : GraphicsViewInput(parent), count(0), scene_(s)
{
  qline = new QGraphicsLineItem();
  qcarc = new CircularArcGraphicsItem<K>();
  qline->hide();
  qcarc->hide();
  s->addItem(qline);
  s->addItem(qcarc);
}


template <typename K>
GraphicsViewCircularArcInput<K>::~GraphicsViewCircularArcInput()
{
  // do not delete qline and qcarc, because 's' owns them and will delete them.
}


template <typename K>
void 
GraphicsViewCircularArcInput<K>::mousePressEvent(QGraphicsSceneMouseEvent *event)
{  
  if(event->modifiers()  & ::Qt::ShiftModifier){
    return;
  }
  if(count == 0){
    qp = event->scenePos();
    p = convert(qp);
    qline->setLine(QLineF(qp, qp));
    qline->show();
    count = 1;
  } else if(count == 1){
    qq = event->scenePos();
    qline->setLine(QLineF(qp, qq));
    q = convert(qq);
    if( (event->button() == ::Qt::RightButton) && (p != q) ){
      qline->hide();
      emit generate(CGAL::make_object(Line_arc_2(Segment_2(p,q))));
      count = 0;
    } else {
      count = 2;
    }
  } else if(count == 2){
    qr  = event->scenePos();
    r = convert(qr);
    typename K::Collinear_2 collinear;
    if(! collinear(p,q,r)){
      qcarc->hide();
      emit generate(make_object(qcarc->arc()));
      count = 0;
    }
  }
}


template <typename K>
void 
GraphicsViewCircularArcInput<K>::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  if(count == 0){
    qcarc->hide();
    qline->hide();
    return;
  } else if(count == 1) {
    qcarc->hide();
    qq = event->scenePos();
    q = convert(qq);
    qline->setLine(QLineF(qp, qq));
    qline->show();
  } else if(count == 2){
    qline->hide();
    qr = event->scenePos();
    r = convert(qr);
    typename K::Collinear_2 collinear;
    if(collinear(p,q,r)){
      qcarc->hide();
      return;
    } else {
      if(CGAL::orientation(p, q, r) == CGAL::RIGHT_TURN) {
	ap = p; ar = r; aq = q;
	qcarc->setArc(Circular_arc_2(p,r,q));
      } else {
	ap = q; ar = r; aq = p;
	qcarc->setArc(Circular_arc_2(q,r,p));
      }
      qcarc->show();
    }
  }
}


template <typename K>
void 
GraphicsViewCircularArcInput<K>::keyPressEvent ( QKeyEvent * event ) 
{
  if(event->key() == ::Qt::Key_Delete){
    if(count>0){
      --count;
    }
  }
  
  if(event->key() == ::Qt::Key_Escape){
    count = 0;
  }
}



template <typename K>
bool 
GraphicsViewCircularArcInput<K>::eventFilter(QObject *obj, QEvent *event)
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

#endif // CGAL_QT_GRAPHICS_VIEW_CIRCULAR_ARC_INPUT_H
