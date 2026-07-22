// Copyright (c) 2024  GeometryFactory Sarl (France).
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

#ifndef CGAL_QT_MULTIPOLYGON_WITH_HOLES_GRAPHICS_ITEM_H
#define CGAL_QT_MULTIPOLYGON_WITH_HOLES_GRAPHICS_ITEM_H

#include <CGAL/license/GraphicsView.h>


#include <CGAL/Bbox_2.h>
#include <CGAL/Multipolygon_with_holes_2.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>
#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename P>
class MultipolygonWithHolesGraphicsItem : public GraphicsItem
{
  typedef typename P::Traits Traits;
public:
  MultipolygonWithHolesGraphicsItem(P* p_);

  void modelChanged();

public:
  QRectF boundingRect() const;

  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);


  const QBrush& brush() const
  {
    return brush_;
  }


  void setBrush(const QBrush& b)
  {
    brush_ = b;
  }

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


protected:
  void updateBoundingBox();

  P * poly;
  QPainter* m_painter;
  PainterOstream<Traits> painterostream;

  typename Traits::Point_2 p;
  QRectF bounding_rect;

  QBrush brush_;
  QPen vertices_pen;
  QPen edges_pen;
  bool draw_vertices;
};


template <typename P>
MultipolygonWithHolesGraphicsItem<P>::MultipolygonWithHolesGraphicsItem(P * p_)
  :  poly(p_), painterostream(0),
     draw_vertices(true)
{
  setVerticesPen(QPen(::Qt::red, 3.));
  if(poly->number_of_polygons_with_holes() == 0){
    this->hide();
  }
  updateBoundingBox();
  setZValue(0);
}

template <typename P>
QRectF
MultipolygonWithHolesGraphicsItem<P>::boundingRect() const
{
  return bounding_rect;
}




template <typename P>
void
MultipolygonWithHolesGraphicsItem<P>::paint(QPainter *painter,
                                            const QStyleOptionGraphicsItem * /*option*/,
                                            QWidget * /*widget*/)
{
  using Polygon_with_holes_2 = typename P::Polygon_with_holes_2;
  using General_polygon_2 = typename Polygon_with_holes_2::General_polygon_2;
  Converter<Traits> convert;

  for(const auto& pwh : poly->polygons_with_holes()){
    QPainterPath border;
    General_polygon_2 boundary =  pwh.outer_boundary();
    painter->setPen(this->edgesPen());

    typename General_polygon_2::Vertex_iterator it = pwh.outer_boundary().vertices_begin();
    QPointF firstPoint = convert(*it);
    border.moveTo(firstPoint);
    painterostream = PainterOstream<Traits>(painter);

    for(++it;
        it != pwh.outer_boundary().vertices_end();
        ++it){
      border.lineTo(convert(*it));
    }
    border.lineTo(firstPoint);


    for(typename Polygon_with_holes_2::Hole_const_iterator hit = pwh.holes_begin();
        hit != pwh.holes_end();
        ++hit){
      typename General_polygon_2::Vertex_iterator it = hit->vertices_begin();
      QPointF firstPoint = convert(*it);
      border.moveTo(firstPoint);
      for(++it;
          it != hit->vertices_end();
          ++it){
        border.lineTo(convert(*it));
      }
      border.lineTo(firstPoint);
    }

    painter->setBrush(this->brush());
    painter->drawPath(border);
  }


  if(drawVertices()) {
    painter->setPen(verticesPen());
    QTransform matrix = painter->worldTransform();
    painter->resetTransform();

    for(const auto& pwh : poly->polygons_with_holes()){
      for(typename General_polygon_2::Vertex_iterator it = pwh.outer_boundary().vertices_begin();
          it != pwh.outer_boundary().vertices_end();
          it++){
        QPointF point = matrix.map(convert(*it));
        painter->drawPoint(point);
      }
      for(typename Polygon_with_holes_2::Hole_const_iterator hit = pwh.holes_begin();
          hit != pwh.holes_end();
          ++hit){
        for(typename General_polygon_2::Vertex_iterator it = hit->vertices_begin();
          it != hit->vertices_end();
          it++){
        QPointF point = matrix.map(convert(*it));
        painter->drawPoint(point);
      }
      }
    }
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename P>
void
MultipolygonWithHolesGraphicsItem<P>::updateBoundingBox()
{
  Converter<Traits> convert;
  prepareGeometryChange();
  if(poly->number_of_polygons_with_holes() == 0){
    return;
  }
  bounding_rect = convert(poly->bbox());
}


template <typename P>
void
MultipolygonWithHolesGraphicsItem<P>::modelChanged()
{
  if((poly->number_of_polygons_with_holes() == 0) ){
    this->hide();
  } else if((poly->number_of_polygons_with_holes() > 0) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_MYLTIPOLYGON_WITH_HOLES_GRAPHICS_ITEM_H
