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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/next/GraphicsView/include/CGAL/Qt/TriangulationGraphicsItem.h $
// $Id: TriangulationGraphicsItem.h 67117 2012-01-13 18:14:48Z lrineau $
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_TRIANGULATION_GRAPHICS_ITEM_H
#define CGAL_QT_TRIANGULATION_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/apply_to_range.h>
#include "HyperbolicPainterOstream.h"
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename T>
class TriangulationGraphicsItem : public GraphicsItem
{
  typedef typename T::Geom_traits Geom_traits;
public:
  TriangulationGraphicsItem(T* t_);

  void modelChanged();

public:

  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  virtual void operator()(typename T::Face_handle fh);

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
  bool visible_edges;
  bool visible_vertices;
};


template <typename T>
TriangulationGraphicsItem<T>::TriangulationGraphicsItem(T * t_)
  :  t(t_), painterostream(0),
     bb(0,0,0,0), bb_initialized(false),
     visible_edges(true), visible_vertices(true)
{
  if(t->number_of_vertices() == 0){
    this->hide();
  }
  updateBoundingBox();
  setZValue(3);
}

template <typename T>
QRectF 
TriangulationGraphicsItem<T>::boundingRect() const
{
  return bounding_rect;
}


template <typename T>
void 
TriangulationGraphicsItem<T>::operator()(typename T::Face_handle fh)
{
  if(visible_edges) {
    for(int i=0; i<3; i++) {
      if(fh < fh->neighbor(i) || !t->is_Delaunay_hyperbolic(fh->neighbor(i))){
        painterostream << t->hyperbolic_segment(fh,i);
      }
    }
  }
  if(visible_vertices) {
    for(int i=0; i<3; i++) {
      paintVertex(fh->vertex(i));
    }
  }
}

template <typename T>
void 
TriangulationGraphicsItem<T>::drawAll(QPainter *painter)
{
  QPen pen;
  pen.setWidthF(0.005);
  pen.setBrush(::Qt::black);
  painter->setPen(edges_pen);
  painterostream = PainterOstream<Geom_traits>(painter);
 
  if(visibleEdges()) {
    for(typename T::All_edges_iterator eit = t->all_edges_begin();
        eit != t->all_edges_end();
        ++eit){
      painterostream << t->hyperbolic_segment(*eit);
    }
  }
  paintVertices(painter);
}

template <typename T>
void 
TriangulationGraphicsItem<T>::paintVertices(QPainter *painter)
{
  if(visibleVertices()) {
    Converter<Geom_traits> convert;

    painter->setPen(vertices_pen);
    QMatrix matrix = painter->matrix();
    painter->resetMatrix();
    for(typename T::All_vertices_iterator it = t->all_vertices_begin();
        it != t->all_vertices_end();
        it++){
      QPointF point = matrix.map(convert(t->point(it)));
      painter->drawPoint(point);
    }
  }
}

template <typename T>
void 
TriangulationGraphicsItem<T>::paintVertex(typename T::Vertex_handle vh)
{
  Converter<Geom_traits> convert;
  m_painter->setPen(vertices_pen);
  QMatrix matrix = m_painter->matrix();
  m_painter->resetMatrix();
  m_painter->drawPoint(matrix.map(convert(t->point(vh))));
  m_painter->setMatrix(matrix);
}

template <typename T>
void 
TriangulationGraphicsItem<T>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * /*widget*/)
{
  if(t->dimension()<2 || option->exposedRect.contains(boundingRect())) {
    drawAll(painter);
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
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename T>
void 
TriangulationGraphicsItem<T>::updateBoundingBox()
{
  prepareGeometryChange();
  bounding_rect = QRectF(-1,-1,2,2);
}


template <typename T>
void 
TriangulationGraphicsItem<T>::modelChanged()
{
  if((t->number_of_vertices() == 0)){
    this->hide();
  } else if((t->number_of_vertices() > 0) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_TRIANGULATION_GRAPHICS_ITEM_H
