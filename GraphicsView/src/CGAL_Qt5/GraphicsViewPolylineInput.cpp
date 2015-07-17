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

#include <QGraphicsItem>
#include <QGraphicsPathItem> 
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QPolygonF>
#include <QPainterPath>
#include <QEvent>
#include <QKeyEvent>

#include <CGAL/Qt/GraphicsViewPolylineInput.h>

namespace CGAL {
namespace Qt {

GraphicsViewPolylineInput_non_templated_base::
GraphicsViewPolylineInput_non_templated_base(QObject* parent,
                                   QGraphicsScene* s,
                                   int n,
                                   bool closed)
  : GraphicsViewInput(parent), closed_(closed), path_item(NULL), b(NULL), e(NULL), n_(n), scene_(s)
{}


bool
GraphicsViewPolylineInput_non_templated_base::mousePressEvent(QGraphicsSceneMouseEvent *event)
{ 
  if( event->modifiers() ){
    return false;
  }
  if( event->button() != ::Qt::RightButton
      && event->button() != ::Qt::LeftButton ){
    return false;
  }
  polygon.push_back(event->scenePos());
  if(path_item){
    scene_->removeItem(path_item);
    delete path_item;
    path_item = NULL;
  }
  if( (event->button() == ::Qt::RightButton) || (polygon.size() == n_) ){
    // call the virtual function generate_polygon(), that emit a
    // CGAL::Object containing a list of points
    generate_polygon();
    polygon.clear();
    if(b){
      scene_->removeItem(b);
      delete b;
      b = NULL;
    }
    if(e){
      scene_->removeItem(e);
      delete e;
      e = NULL;
    }
    return true;
  }
  if(event->button() == ::Qt::LeftButton){
    QPainterPath qpp;
    qpp.addPolygon(polygon);
    path_item = new QGraphicsPathItem(qpp);
    path_item->setPen(QPen(::Qt::red, 0, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
    scene_->addItem(path_item);
    return true;
  }
  return false;
}


void 
GraphicsViewPolylineInput_non_templated_base::rubberbands(const QPointF& p)
{
  if(polygon.empty()){
    return;
  }
  if(!b && closed_ ){
    b = new QGraphicsLineItem();
    b->setPen(QPen(::Qt::red, 0, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
    scene_->addItem(b);
  }
  if( !e){
    e = new QGraphicsLineItem();    
    e->setPen(QPen(::Qt::red, 0, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
    scene_->addItem(e);
  }
  if(closed_){
    QLineF bLine(polygon.front(), p);
    b->setLine(bLine);
  }
  QLineF eLine(polygon.back(), p);
  e->setLine(eLine); 
}


void 
GraphicsViewPolylineInput_non_templated_base::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
  sp = event->scenePos();
  rubberbands(sp);
}


bool
GraphicsViewPolylineInput_non_templated_base::keyPressEvent ( QKeyEvent * event ) 
{
  if( event->modifiers() )
    return false;

  switch(event->key())
  {
  case ::Qt::Key_Delete:
  case ::Qt::Key_Escape:
  case ::Qt::Key_Backspace:
    break;
  default:
    return false;
  }
  if(polygon.empty()){
    return true;
  }
  polygon.pop_back();
  if(polygon.empty()){
    if(b){
      scene_->removeItem(b);
      delete b;
      b = NULL;
    }
    if(e){
      scene_->removeItem(e);
      delete e;
      e = NULL;
    }
    return true;
  }
  if(path_item){
    scene_->removeItem(path_item);
    delete path_item;
    path_item = NULL;
  }
  QPainterPath qpp;
  qpp.addPolygon(polygon);
  path_item = new QGraphicsPathItem(qpp);
  path_item->setPen(QPen(::Qt::red, 0, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
  scene_->addItem(path_item);
  rubberbands(sp);
  return true;
}



bool 
GraphicsViewPolylineInput_non_templated_base::eventFilter(QObject *obj, QEvent *event)
{
  if (event->type() == QEvent::GraphicsSceneMousePress) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    if(!mousePressEvent(mouseEvent)) {
      // standard event processing if mousePressEvent has returned false
      return QObject::eventFilter(obj, event);
    }
  } else if (event->type() == QEvent::GraphicsSceneMouseMove) {
    QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
    mouseMoveEvent(mouseEvent);
    return QObject::eventFilter(obj, event);
  } else if (event->type() == QEvent::KeyPress) {
    QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);
    if(!keyPressEvent(keyEvent)) {
      return QObject::eventFilter(obj, event);
    }
  }
  // standard event processing if keyPressEvent has returned false
  return QObject::eventFilter(obj, event);
}

} // namespace Qt
} // namespace CGAL
