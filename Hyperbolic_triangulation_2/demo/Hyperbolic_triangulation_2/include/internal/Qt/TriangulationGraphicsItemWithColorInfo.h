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
  setVerticesPen(QPen(::Qt::red, 3.));
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
      if(fh < fh->neighbor(i) || t->is_infinite(fh->neighbor(i))){
        m_painter->setPen(this->edgesPen());
        painterostream << t->segment(fh,i);
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
  //delete
  QPen temp = painter->pen();
  QPen old = temp;
  temp.setWidthF(/*0.0035*/0.0045);
  painter->setPen(temp);
  //

  painterostream = PainterOstream<Geom_traits>(painter);

  if(visibleEdges()) {
    for(typename T::Finite_edges_iterator eit = t->finite_edges_begin();
        eit != t->finite_edges_end();
        ++eit){
      painterostream << t->segment(*eit);
    }
  }

  //delete
  painter->setPen(old);
  //

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

      // draw vertices with color storing in their info
      if(it->info().getColor() == 0) {
        painter->setPen(QPen(::Qt::red, 3.));
      }

      if(it->info().getColor() == 1) {
        painter->setPen(QPen(::Qt::green, 3.));
      }

      if(it->info().getColor() == 2) {
        painter->setPen(QPen(::Qt::cyan, 3.));
      }

      if(it->info().getColor() == 3) {
        painter->setPen(QPen(::Qt::magenta, 3.));
      }

      if(it->info().getColor() == 6) {
        painter->setPen(QPen(::Qt::yellow, 3.));
      }

      if(it->info().getColor() == 5) {
        // brown
        painter->setPen(QPen(QColor(139, 69, 19), 3.));
      }

      if(it->info().getColor() == 4) {
        painter->setPen(QPen(::Qt::blue, 3.));
      }


      if(it->info().getColor() == 7) {
        // orange
        QColor orange = QColor(255, 165, 0);
        painter->setPen(QPen(orange, 3.));
      }

      if(it->info().getColor() == 8) {
        // dark green
        QColor blue = QColor(0, 102, 51);
        painter->setPen(QPen(blue, 3.));
      }

      if(it->info().getColor() == 9) {
        // purple
        QColor blue = QColor(102, 0, 102);
        painter->setPen(QPen(blue, 3.));
      }

      if(it->info().getColor() == 10) {
        // close to blue
        QColor blue = QColor(131, 111, 255);
        painter->setPen(QPen(blue, 3.));
      }

      //

      // delete
      QPen temp = painter->pen();
      QPen old = temp;
      temp.setWidth(9);

      double px = to_double(t->point(it).x());
      double py = to_double(t->point(it).y());
      double dist = px*px + py*py;
      if(dist > 0.25) {
        temp.setWidth(8);//6
      }
      if(dist > 0.64) {
        temp.setWidth(7);//5
      }
      if(dist > 0.81) {
        temp.setWidth(5);//3
      }
      if(dist > 0.92) {
        temp.setWidth(4);//3
      }
      if(dist > 0.98) {
        temp.setWidth(3);//3
      }
      //painter->setPen(temp);

      QPointF point = matrix.map(convert(t->point(it)));
      painter->drawPoint(point);

      painter->setPen(old);

      /*
      QBrush temp = painter->brush();
      QBrush old = temp;
      temp.setColor(painter->pen().color());

      painter->setBrush(temp);

      painter->drawEllipse(point, 10, 10);

      painter->setBrush(old);
       */
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
  m_painter->drawPoint(matrix.map(convert(t->point(vh))));
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
  if(t->number_of_vertices() == 0){
    bb = Bbox_2(0,0,0,0);
    bb_initialized = false;
    return;
  } else if(! bb_initialized){
    bb = t->point(t->finite_vertices_begin()).bbox();
    bb_initialized = true;
  }

  if(t->dimension() <2){
    for(typename T::Finite_vertices_iterator it = t->finite_vertices_begin();
        it != t->finite_vertices_end();
        ++it){
      bb = bb + t->point(it).bbox();
    }
  } else {
    typename T::Vertex_handle inf = t->infinite_vertex();
    typename T::Vertex_circulator vc = t->incident_vertices(inf), done(vc);
    do {
      bb = bb + t->point(vc).bbox();
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
