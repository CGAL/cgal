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

#ifndef CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
#define CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
//#include <CGAL/apply_to_range.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename T>
class ArrangementGraphicsItem : public GraphicsItem
{
  typedef typename T::Traits_2 Geom_traits;

public:
  ArrangementGraphicsItem(T* t_);

  void modelChanged();

public:

  QRectF boundingRect() const;
  
  virtual void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  //virtual void operator()(typename T::Face_handle fh);

  const QPen& verticesPen() const
  {
    return this->vertices_pen;
  }

  const QPen& edgesPen() const
  {
    return this->edges_pen;
  }

  void setVerticesPen(const QPen& pen)
  {
    this->vertices_pen = pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    this->edges_pen = pen;
  }

  bool visibleVertices() const
  {
    return this->visible_vertices;
  }

  void setVisibleVertices(const bool b)
  {
    this->visible_vertices = b;
    this->update();
  }

  bool visibleEdges() const
  {
    return this->visible_edges;
  }

  void setVisibleEdges(const bool b)
  {
    this->visible_edges = b;
    this->update();
  }

protected:
  //virtual void drawAll(QPainter *painter);
  //void paintVertices(QPainter *painter);
  //void paintOneVertex(const typename T::Point& point);
  //virtual void paintVertex(typename T::Vertex_handle vh);
  void updateBoundingBox();

  T * t;
  QPainter* m_painter;
  PainterOstream<Geom_traits> painterostream;

  typename T::Vertex_handle vh;
  typename Geom_traits::Point_2 p;
  CGAL::Bbox_2 bb;  
  bool bb_initialized;
  QRectF bounding_rect;

  QPen vertices_pen;
  QPen edges_pen;
  bool visible_edges;
  bool visible_vertices;
};


template <typename T>
ArrangementGraphicsItem<T>::ArrangementGraphicsItem(T * t_)
  :  t(t_), painterostream(0),
     bb(0,0,0,0), bb_initialized(false),
     visible_edges(true), visible_vertices(true)
{
  this->setVerticesPen(QPen(::Qt::red, 3.));
  if (t->number_of_vertices() == 0) {
    this->hide();
  }
  this->updateBoundingBox();
  this->setZValue(3);
}

template <typename T>
QRectF 
ArrangementGraphicsItem<T>::boundingRect() const
{
  return this->bounding_rect;
}

/*
template <typename T>
void 
TriangulationGraphicsItem<T>::operator()(typename T::Face_handle fh)
{
  if (this->visible_edges) {
    for (int i=0; i<3; i++) {
      if (fh < fh->neighbor(i) || this->t->is_infinite(fh->neighbor(i))) {
        this->m_painter->setPen(this->edgesPen());
        this->painterostream << this->t->segment(fh, i);
      }
    }
  }
  if(visible_vertices) {
    for (int i=0; i<3; i++) {
      this->paintVertex(fh->vertex(i));
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
TriangulationGraphicsItem<T>::paintOneVertex(const typename T::Point& point)
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
TriangulationGraphicsItem<T>::paintVertex(typename T::Vertex_handle vh)
{
  Converter<Geom_traits> convert;

  m_painter->setPen(this->verticesPen());
  QMatrix matrix = m_painter->matrix();
  m_painter->resetMatrix();
  m_painter->drawPoint(matrix.map(convert(vh->point())));
  m_painter->setMatrix(matrix);
}
*/

template <typename T>
void 
ArrangementGraphicsItem<T>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * /*widget*/)
{
    std::cout << "TriangulationGraphicsItem::paint not yet implemented" << std::endl;
#if 0
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
#endif
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename T>
void 
ArrangementGraphicsItem<T>::updateBoundingBox()
{
    std::cout << "updateBoundingBox stub" << std::endl;
#if 0
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
#endif
}


template <typename T>
void 
ArrangementGraphicsItem<T>::modelChanged()
{
  /*
  if((this->t->is_empty()) ){
    this->hide();
  } else if((t->number_of_vertices() > 0) && (! this->isVisible())){
    this->show();
  }
  */
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
