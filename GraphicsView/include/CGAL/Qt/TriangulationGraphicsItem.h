// Copyright (c) 2008  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_TRIANGULATION_GRAPHICS_ITEM_H
#define CGAL_QT_TRIANGULATION_GRAPHICS_ITEM_H

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

  const QPen& verticesPen() const
  {
    return vertices_pen;
  }

  QPen& verticesPen()
  {
    return vertices_pen;
  }
  const QPen& edgesPen() const
  {
    return edges_pen;
  }

  QPen& edgesPen()
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
  setVerticesPen(QPen(::Qt::red, 4.));
  setEdgesPen(QPen(::Qt::black, 0, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
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
TriangulationGraphicsItem<T>::drawAll(QPainter *painter)
{
  painterostream = PainterOstream<Geom_traits>(painter);

  if(visibleEdges()) {
    for(typename T::Finite_edges_iterator eit = t->finite_edges_begin();
        eit != t->finite_edges_end();
        ++eit){
      painterostream << t->segment(*eit);
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

    painter->setPen(verticesPen());
    QTransform matrix = painter->worldTransform();
    painter->resetTransform();
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
TriangulationGraphicsItem<T>::paintOneVertex(const typename T::Point& point)
{
  Converter<Geom_traits> convert;

  m_painter->setPen(this->verticesPen());
  QTransform matrix = m_painter->worldTransform();
  m_painter->resetTransform();
  m_painter->drawPoint(matrix.map(convert(point)));
  m_painter->setWorldTransform(matrix);
}

template <typename T>
void
TriangulationGraphicsItem<T>::paintVertex(typename T::Vertex_handle vh)
{
  Converter<Geom_traits> convert;

  m_painter->setPen(this->verticesPen());
  QTransform matrix = m_painter->worldTransform();
  m_painter->resetTransform();
  m_painter->drawPoint(matrix.map(convert(vh->point())));
  m_painter->setWorldTransform(matrix);
}

template <typename T>
void
TriangulationGraphicsItem<T>::paint(QPainter *painter,
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * /*widget*/)
{
  painter->setPen(this->edgesPen());
//   painter->drawRect(boundingRect());
  if ( t->dimension()<2 || option->exposedRect.contains(boundingRect()) ) {
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
TriangulationGraphicsItem<T>::modelChanged()
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

#endif // CGAL_QT_TRIANGULATION_GRAPHICS_ITEM_H
