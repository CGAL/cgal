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

#ifndef CGAL_QT_REGULAR_TRIANGULATION_GRAPHICS_ITEM_H
#define CGAL_QT_REGULAR_TRIANGULATION_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/apply_to_range.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename T>
class RegularTriangulationGraphicsItem : public GraphicsItem
{
  typedef typename T::Geom_traits Geom_traits;      
  typedef typename Kernel_traits<typename T::Bare_point>::Kernel K;
  typedef typename K::Segment_2 Segment_2;

public:
  RegularTriangulationGraphicsItem(T* t_);

  void modelChanged();

public:

  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  virtual void operator()(typename T::Face_handle fh);

  const QPen& verticesPen() const
  {
    return vertices_pen;
  }

  const QPen& edgesPen() const
  {
    return edges_pen;
  }

  void setVerticesPen(const QPen& pen)
  {
    vertices_pen = pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    edges_pen = pen;
  }

  bool visibleVertices() const
  {
    return visible_vertices;
  }

  void setVisibleVertices(const bool b)
  {
    visible_vertices = b;
    update();
  }

  bool visibleEdges() const
  {
    return visible_edges;
  }

  void setVisibleEdges(const bool b)
  {
    visible_edges = b;
    update();
  }

protected:
  virtual void drawAll(QPainter *painter);
  void paintVertices(QPainter *painter);
  void paintOneVertex(const typename T::Point& point);
  void updateBoundingBox();

  T * t;
  QPainter* m_painter;

  PainterOstream<K> painterostream;

  typename T::Vertex_handle vh;
  typename T::Point p;
  CGAL::Bbox_2 bb;  
  bool bb_initialized;
  QRectF bounding_rect;

  QPen vertices_pen;
  QPen edges_pen;
  bool visible_edges;
  bool visible_vertices;
};


template <typename T>
RegularTriangulationGraphicsItem<T>::RegularTriangulationGraphicsItem(T * t_)
  :  t(t_), painterostream(0), bb(0,0,0,0), bb_initialized(false),
     visible_edges(true), visible_vertices(true)
     
{
  setVerticesPen(QPen(::Qt::red, 3.));
  if(t->number_of_vertices() == 0){
    this->hide();
  }
  updateBoundingBox();
  setZValue(3);
}

template <typename T>
QRectF 
RegularTriangulationGraphicsItem<T>::boundingRect() const
{
  return bounding_rect;
}


template <typename T>
void 
RegularTriangulationGraphicsItem<T>::operator()(typename T::Face_handle fh)
{
  
  if(visible_edges) {
    for (int i=0; i<3; i++) {
      if (fh < fh->neighbor(i) || t->is_infinite(fh->neighbor(i))){
        m_painter->setPen(this->edgesPen());
	Segment_2 s(fh->vertex(T::cw(i))->point().point(), fh->vertex(T::ccw(i))->point().point());
        painterostream << s;
      }
    }
  }
  if(visible_vertices) {
    for (int i=0; i<3; i++) {
      paintOneVertex(fh->vertex(i)->point());
    }
  }
  
}

template <typename T>
void 
RegularTriangulationGraphicsItem<T>::drawAll(QPainter *painter)
{
  painterostream = PainterOstream<K>(painter);
 
  
  if(visibleEdges()) {
    for(typename T::Finite_edges_iterator eit = t->finite_edges_begin();
        eit != t->finite_edges_end();
        ++eit){
      typename T::Face_handle fh = eit->first;
      int i = eit->second;
      Segment_2 s(fh->vertex(T::cw(i))->point().point(), fh->vertex(T::ccw(i))->point().point());
      painterostream << s;
    }
  }
  paintVertices(painter);
  
}

template <typename T>
void 
RegularTriangulationGraphicsItem<T>::paintVertices(QPainter *painter)
{
  if(visibleVertices()) {
      Converter<K> convert;

    painterostream = PainterOstream<K>(painter);
    painter->setPen(verticesPen());
    
    for(typename T::Finite_vertices_iterator it = t->finite_vertices_begin();
        it != t->finite_vertices_end();
        it++){

      typename K::Circle_2 circ(it->point().point(), it->point().weight()); 
      painterostream << circ;
    }


    QMatrix matrix = painter->matrix();
    painter->resetMatrix();
    for(typename T::Finite_vertices_iterator it = t->finite_vertices_begin();
        it != t->finite_vertices_end();
        it++){
      QPointF point = matrix.map(convert(it->point().point()));
      painter->drawPoint(point);
    }
  }
}

template <typename T>
void 
RegularTriangulationGraphicsItem<T>::paintOneVertex(const typename T::Point& point)
{
  Converter<K> convert;

  m_painter->setPen(this->verticesPen());
  QMatrix matrix = m_painter->matrix();
  m_painter->resetMatrix();
  m_painter->drawPoint(matrix.map(convert(point)));
  m_painter->setMatrix(matrix);
}

template <typename T>
void 
RegularTriangulationGraphicsItem<T>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem * /*option*/,
                                    QWidget * /*widget*/)
{
  painter->setPen(this->edgesPen());
  drawAll(painter);
  return;
#if 0
    if ( t->dimension()<2 || option->exposedRect.contains(boundingRect()) ) {
    drawAll(painter);
  } else {
    m_painter = painter;
    painterostream = PainterOstream<K>(painter);
    CGAL::apply_to_range (*t, 
                          typename T::Point(typename T::Bare_point(option->exposedRect.left(),
								   option->exposedRect.bottom()),0), 
                          typename T::Point(typename T::Bare_point(option->exposedRect.right(),
								   option->exposedRect.top()),0), 
                          *this);
  }
#endif
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename T>
void 
RegularTriangulationGraphicsItem<T>::updateBoundingBox()
{
  prepareGeometryChange();
  if(t->number_of_vertices() == 0){
    bb = Bbox_2(0,0,0,0);
    bb_initialized = false;
    return;
  } else if(! bb_initialized){
    bb = t->finite_vertices_begin()->point().bbox();
    bb_initialized = true;
  }
  
  if(t->dimension() <2){
    for(typename T::Finite_vertices_iterator it = t->finite_vertices_begin();
	it != t->finite_vertices_end();
	++it){
      typename K::Circle_2 circ(it->point().point(), it->point().weight());
      bb = bb + circ.bbox();
    }
  } else {
    typename T::Vertex_handle inf = t->infinite_vertex();
    typename T::Vertex_circulator vc = t->incident_vertices(inf), done(vc);
    do {
      typename K::Circle_2 circ(vc->point().point(), vc->point().weight());
      bb = bb + circ.bbox();
      ++vc;
    } while(vc != done);
  }
  bounding_rect = QRectF(bb.xmin(),
                         bb.ymin(),
                         bb.xmax()-bb.xmin(),
                         bb.ymax()-bb.ymin());
}


template <typename T>
void 
RegularTriangulationGraphicsItem<T>::modelChanged()
{
  if((t->number_of_vertices() == 0) ){
    this->hide();
  } else if((t->number_of_vertices() > 0) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_REGULAR_TRIANGULATION_GRAPHICS_ITEM_H
