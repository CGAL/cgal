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

#ifndef CGAL_QT_ALPHA_SHAPE_GRAPHICS_ITEM_H
#define CGAL_QT_ALPHA_SHAPE_GRAPHICS_ITEM_H

#include <CGAL/license/GraphicsView.h>


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
class AlphaShapeGraphicsItem : public GraphicsItem
{
  typedef typename T::Geom_traits Geom_traits;
public:
  AlphaShapeGraphicsItem(T* t_);

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

  const QPen& regularEdgesPen() const
  {
    return regular_edges_pen;
  }

  const QPen& singularEdgesPen() const
  {
    return singular_edges_pen;
  }

  const QBrush& regularFacesBrush() const
  {
    return regular_faces_brush;
  }

  void setVerticesPen(const QPen& pen)
  {
    vertices_pen = pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    edges_pen = pen;
  }

  void setRegularEdgesPen(const QPen& pen)
  {
    regular_edges_pen = pen;
  }

  void setSingularEdgesPen(const QPen& pen)
  {
    singular_edges_pen = pen;
  }

  void setRegularFacesBrush(const QBrush& b)
  {
    regular_faces_brush = b;
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
  virtual void paintVertex(typename T::Vertex_handle vh);
  void updateBoundingBox();

  T * t;
  QPainter* m_painter;
  PainterOstream<Geom_traits> painterostream;

  typename T::Vertex_handle vh;
  typename T::Point p;
  CGAL::Bbox_2 bb;  
  bool bb_initialized;
  QRectF bounding_rect;

  QPen vertices_pen;
  QPen edges_pen;
  QPen regular_edges_pen;
  QPen singular_edges_pen;
  QBrush regular_faces_brush;
  bool visible_edges;
  bool visible_vertices;
};


template <typename T>
AlphaShapeGraphicsItem<T>::AlphaShapeGraphicsItem(T * t_)
  :  t(t_), painterostream(0),
     bb(0,0,0,0), bb_initialized(false),
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
AlphaShapeGraphicsItem<T>::boundingRect() const
{
  return bounding_rect;
}


template <typename T>
void 
AlphaShapeGraphicsItem<T>::operator()(typename T::Face_handle fh)
{
  if(visible_edges) {
    for (int i=0; i<3; i++) {
      if (fh < fh->neighbor(i) || t->is_infinite(fh->neighbor(i))){
        m_painter->setPen(this->edgesPen());
        painterostream << t->segment(fh,i);
      }
    }
  }
  if(visible_vertices) {
    for (int i=0; i<3; i++) {
      paintVertex(fh->vertex(i));
    }
  }
}

template <typename T>
void 
AlphaShapeGraphicsItem<T>::drawAll(QPainter *painter)
{
  painterostream = PainterOstream<Geom_traits>(painter);

  painter->setBrush(regularFacesBrush());
  typedef typename T::Gt::Triangle_2 Triangle_2;
  for(typename T::Finite_faces_iterator fit = t->finite_faces_begin();
      fit != t->finite_faces_end();
      ++fit){
    
    if(t->classify(fit) == T::INTERIOR){
      Triangle_2 triangle = t->triangle(fit);
      painterostream << triangle;
    }
  }

  if(visibleEdges()) {
    for(typename T::Finite_edges_iterator eit = t->finite_edges_begin();
        eit != t->finite_edges_end();
        ++eit){
      switch  (t->classify(*eit)) {
      case T::REGULAR :
        painter->setPen(regularEdgesPen());
        break;
      case T::SINGULAR :
        painter->setPen(singularEdgesPen());
        break;
      default:
        painter->setPen(edgesPen());
      }
      painterostream << t->segment(*eit);
    }
  }
  paintVertices(painter);
}

template <typename T>
void 
AlphaShapeGraphicsItem<T>::paintVertices(QPainter *painter)
{
  if(visibleVertices()) {
    Converter<Geom_traits> convert;

    painter->setPen(verticesPen());
    QMatrix matrix = painter->matrix();
    painter->resetMatrix();
    for(typename T::Finite_vertices_iterator it = t->finite_vertices_begin();
        it != t->finite_vertices_end();
        it++){
      QPointF point = matrix.map(convert(it->point()));
      painter->drawPoint(point);
    }
  }
}

template <typename T>
void 
AlphaShapeGraphicsItem<T>::paintOneVertex(const typename T::Point& point)
{
  Converter<Geom_traits> convert;

  m_painter->setPen(this->verticesPen());
  QMatrix matrix = m_painter->matrix();
  m_painter->resetMatrix();
  m_painter->drawPoint(matrix.map(convert(point)));
  m_painter->setMatrix(matrix);
}

template <typename T>
void 
AlphaShapeGraphicsItem<T>::paintVertex(typename T::Vertex_handle vh)
{
  Converter<Geom_traits> convert;

  m_painter->setPen(this->verticesPen());
  QMatrix matrix = m_painter->matrix();
  m_painter->resetMatrix();
  m_painter->drawPoint(matrix.map(convert(vh->point())));
  m_painter->setMatrix(matrix);
}

template <typename T>
void 
AlphaShapeGraphicsItem<T>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem * /*option*/,
                                    QWidget * /*widget*/)
{
  painter->setPen(this->edgesPen());
//   painter->drawRect(boundingRect());
  // if ( t->dimension()<2 || option->exposedRect.contains(boundingRect()) ) {
    drawAll(painter);
    /*
  } else {
    m_painter = painter;
    painterostream = PainterOstream<Geom_traits>(painter);
    CGAL::apply_to_range (*t, 
                          typename T::Point(option->exposedRect.left(),
                                            option->exposedRect.bottom()), 
                          typename T::Point(option->exposedRect.right(),
                                            option->exposedRect.top()), 
                          *this);
  }
    */
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename T>
void 
AlphaShapeGraphicsItem<T>::updateBoundingBox()
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
      bb = bb + it->point().bbox();
    }
  } else {
    typename T::Vertex_handle inf = t->infinite_vertex();
    typename T::Vertex_circulator vc = t->incident_vertices(inf), done(vc);
    do {
      bb = bb + vc->point().bbox();
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
AlphaShapeGraphicsItem<T>::modelChanged()
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

#endif // CGAL_QT_ALPHA_SHAPE_GRAPHICS_ITEM_H
