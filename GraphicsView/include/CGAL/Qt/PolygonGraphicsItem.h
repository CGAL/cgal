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

#ifndef CGAL_QT_POLYGON_GRAPHICS_ITEM_H
#define CGAL_QT_POLYGON_GRAPHICS_ITEM_H

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

template <typename P>
class PolygonGraphicsItem : public GraphicsItem
{
  typedef typename P::Traits Traits;
public:
  PolygonGraphicsItem(P* p_);

  void modelChanged();

public:
  QRectF boundingRect() const;

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);


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

  bool drawVertices() const
  {
    return draw_vertices;
  }

  void setDrawVertices(const bool b)
  {
    draw_vertices = b;
    update();
  }

  bool drawEdges() const
  {
    return draw_edges;
  }

  void setDrawEdges(const bool b)
  {
    draw_edges = b;
    update();
  }

protected:
  void updateBoundingBox();

  P * poly;
  QPainter* m_painter;
  PainterOstream<Traits> painterostream;

  typename P::Point_2 p;
  QRectF bounding_rect;

  QPen vertices_pen;
  QPen edges_pen;
  bool draw_edges;
  bool draw_vertices;
};


template <typename P>
PolygonGraphicsItem<P>::PolygonGraphicsItem(P * p_)
  :  poly(p_), painterostream(0),
     draw_edges(true), draw_vertices(true)
{
  setVerticesPen(QPen(::Qt::red, 3.));
  setEdgesPen(QPen(::Qt::black, 0));
  updateBoundingBox();
  setZValue(3);
}

template <typename P>
QRectF
PolygonGraphicsItem<P>::boundingRect() const
{
  return bounding_rect;
}




template <typename P>
void
PolygonGraphicsItem<P>::paint(QPainter *painter,
                                    const QStyleOptionGraphicsItem * /*option*/,
                                    QWidget * /*widget*/)
{
  painter->setPen(this->edgesPen());
  painterostream = PainterOstream<Traits>(painter);
  if(drawEdges()) {
    for(typename P::Edge_const_iterator eit = poly->edges_begin();
        eit != poly->edges_end();
        ++eit){
      painterostream << *eit;
    }
  }

  if(drawVertices()) {
    Converter<Traits> convert;

    painter->setPen(verticesPen());
    QTransform matrix = painter->worldTransform();
    painter->resetTransform();
    for(typename P::Vertex_iterator it = poly->vertices_begin();
        it != poly->vertices_end();
        it++){
      QPointF point = matrix.map(convert(*it));
      painter->drawPoint(point);
    }
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename P>
void
PolygonGraphicsItem<P>::updateBoundingBox()
{
  Converter<Traits> convert;
  prepareGeometryChange();
  if(poly->size() == 0){
    return;
  }
  bounding_rect = convert(poly->bbox());
}


template <typename P>
void
PolygonGraphicsItem<P>::modelChanged()
{
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_POLYGON_GRAPHICS_ITEM_H
