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

#ifndef CGAL_QT_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H
#define CGAL_QT_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H

#include <CGAL/Bbox_2.h>
#include <CGAL/Polygon_with_holes_2.h>
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
class PolygonWithHolesGraphicsItem : public GraphicsItem
{
  typedef typename P::General_polygon_2::Traits Traits;
public:
  PolygonWithHolesGraphicsItem(P* p_);

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

  typename P::General_polygon_2::Point_2 p;
  QRectF bounding_rect;

  QBrush brush_;
  QPen vertices_pen;
  QPen edges_pen;
  bool draw_vertices;
};


template <typename P>
PolygonWithHolesGraphicsItem<P>::PolygonWithHolesGraphicsItem(P * p_)
  :  poly(p_), painterostream(0),
     draw_vertices(true)   
{
  setVerticesPen(QPen(::Qt::red, 3.));
  if(poly->outer_boundary().size() == 0){
    this->hide();
  }
  updateBoundingBox();
  setZValue(0);
}

template <typename P>
QRectF 
PolygonWithHolesGraphicsItem<P>::boundingRect() const
{
  return bounding_rect;
}




template <typename P>
void 
PolygonWithHolesGraphicsItem<P>::paint(QPainter *painter, 
				       const QStyleOptionGraphicsItem * /*option*/,
				       QWidget * /*widget*/)
{
  Converter<Traits> convert;
  QPainterPath border;
  typename P::General_polygon_2 boundary =  poly->outer_boundary();
  painter->setPen(this->edgesPen());

  typename P::General_polygon_2::Vertex_iterator it = poly->outer_boundary().vertices_begin();
  QPointF firstPoint = convert(*it);
  border.moveTo(firstPoint);
  painterostream = PainterOstream<Traits>(painter);

  for(++it;
      it != poly->outer_boundary().vertices_end();
      ++it){
    border.lineTo(convert(*it)); 
  }
  border.lineTo(firstPoint);

 
  for(typename P::Hole_const_iterator hit = poly->holes_begin();
      hit != poly->holes_end();
      ++hit){
    typename P::General_polygon_2::Vertex_iterator it = hit->vertices_begin();
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

  if(drawVertices()) {

    painter->setPen(verticesPen());
    QMatrix matrix = painter->matrix();
    painter->resetMatrix();
    for(typename P::General_polygon_2::Vertex_iterator it = poly->outer_boundary().vertices_begin();
        it != poly->outer_boundary().vertices_end();
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
PolygonWithHolesGraphicsItem<P>::updateBoundingBox()
{
  Converter<Traits> convert;
  prepareGeometryChange();
  if(poly->outer_boundary().size() == 0){
    return;
  }
  bounding_rect = convert(poly->outer_boundary().bbox());
}


template <typename P>
void 
PolygonWithHolesGraphicsItem<P>::modelChanged()
{
  if((poly->outer_boundary().size() == 0) ){
    this->hide();
  } else if((poly->outer_boundary().size() > 0) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_POLYGON_WITH_HOLES_GRAPHICS_ITEM_H
