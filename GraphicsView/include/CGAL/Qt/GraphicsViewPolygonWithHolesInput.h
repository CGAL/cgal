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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_GRAPHICS_VIEW_POLYGON_WITH_HOLES_INPUT_H
#define CGAL_QT_GRAPHICS_VIEW_POLYGON_WITH_HOLES_INPUT_H

#include <CGAL/license/GraphicsView.h>


#include <list>

#include <QGraphicsView>
#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QKeyEvent>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Qt/GraphicsViewInput.h>
#include <CGAL/Qt/PolygonWithHolesGraphicsItem.h>
#include <CGAL/Qt/GraphicsViewPolylineInput.h>
#include <CGAL/array.h>

namespace CGAL {
namespace Qt {

  /*
    We store a polygon with holes and display it with a graphics item
    We use a PolygonInput tool for entering the boundary and the holes
    We forward most events directly to the polygon input tool
    We only deal with events when the polygon input is not active
    - left click: enter new polygon
    - right click: return result
    - backspace: delete last polygon
    - esc: return with empty result

    todo: check that polygons don't intersect
   */



template <typename K>
class GraphicsViewPolygonWithHolesInput : public GraphicsViewInput
{
public:
  GraphicsViewPolygonWithHolesInput(QObject *parent, QGraphicsScene* s); 
  ~GraphicsViewPolygonWithHolesInput();
  
public Q_SLOTS:
  void processInput(CGAL::Object o);

typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes;

protected:
    
  virtual void keyPressEvent(QKeyEvent *event);
  
  bool eventFilter(QObject *obj, QEvent *event);
  
private:

  Polygon polygon;
  std::list<Polygon> holes; 
  Polygon_with_holes pwh;  // this one collects the input polygons

  CGAL::Qt::PolygonWithHolesGraphicsItem<Polygon_with_holes> * pwhItem;
  CGAL::Qt::GraphicsViewPolylineInput<K> * pi;

  bool polygon_input;
  typedef typename K::Point_2 Point_2;
  QGraphicsScene *scene_;  
};


template <typename K>
GraphicsViewPolygonWithHolesInput<K>::GraphicsViewPolygonWithHolesInput(QObject *parent, QGraphicsScene* s)
  : GraphicsViewInput(parent), scene_(s), polygon_input(false)
{
  pwhItem = new CGAL::Qt::PolygonWithHolesGraphicsItem<Polygon_with_holes>(&pwh);
  pwhItem->setBrush(::Qt::yellow);
  scene_->addItem(pwhItem);
  pwhItem->hide();
  
  pi = new CGAL::Qt::GraphicsViewPolylineInput<K>(parent,s);
  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));

  QObject::connect(this, SIGNAL(modelChanged()),
		   pwhItem, SLOT(modelChanged()));

}

template <typename K>
GraphicsViewPolygonWithHolesInput<K>::~GraphicsViewPolygonWithHolesInput()
{
  //delete pwhItem;
  //delete pi;
}


template <typename K>
void
GraphicsViewPolygonWithHolesInput<K>::processInput(CGAL::Object o)
{
   std::list<Point_2> points;
  if(CGAL::assign(points, o)){
    if((points.size() == 1)&& polygon.size()>0){
    
    } else {
      polygon.clear();
      if(points.front() == points.back()){
	points.pop_back();
      }
      polygon.insert(polygon.vertices_begin(), points.begin(), points.end());
      if(holes.empty()){
	if(polygon.orientation() == CGAL::CLOCKWISE){
	  polygon.reverse_orientation();
	}
      } else {
	if(polygon.orientation() == CGAL::COUNTERCLOCKWISE){
	  polygon.reverse_orientation();
	}
      }
      holes.push_back(polygon);
      typename std::list<Polygon>::iterator it = holes.begin();
      it++;
      pwh = Polygon_with_holes(holes.front(), it, holes.end());
    }
    Q_EMIT( modelChanged());
    polygon_input = false;
  } 
}


template <typename K>
void 
GraphicsViewPolygonWithHolesInput<K>::keyPressEvent ( QKeyEvent * event ) 
{
}



template <typename K>
bool 
GraphicsViewPolygonWithHolesInput<K>::eventFilter(QObject *obj, QEvent *event)
{
  if(polygon_input){
    return pi->eventFilter(obj, event);
  } else {
    if (event->type() == QEvent::GraphicsSceneMousePress) {
      QGraphicsSceneMouseEvent *mouseEvent = static_cast<QGraphicsSceneMouseEvent *>(event);
      
      if(mouseEvent->modifiers()  & ::Qt::ShiftModifier){
	return QObject::eventFilter(obj, event);;
      }
      if(mouseEvent->button() == ::Qt::LeftButton) {
	polygon_input = true;
	return pi->eventFilter(obj, event);
      } else if(mouseEvent->button() == ::Qt::RightButton) {
	Q_EMIT( generate(CGAL::make_object(pwh)));
	pwh.clear();
	holes.clear();
	polygon_input = false;
	Q_EMIT( modelChanged());
      }
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
} 

} // namespace Qt

} // namespace CGAL

#endif // CGAL_QT_GRAPHICS_VIEW_POLYGON_WITH_HOLES_INPUT_H
