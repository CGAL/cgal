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

#ifndef CGAL_QT_POLYLINES_GRAPHICS_ITEM_H
#define CGAL_QT_POLYLINES_GRAPHICS_ITEM_H

#include <CGAL/license/GraphicsView.h>


#include <CGAL/Bbox_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename P>
class PolylinesGraphicsItem : public GraphicsItem
{
  typedef typename P::value_type Polyline;
  typedef typename Polyline::value_type Point_2;
  typedef typename CGAL::Kernel_traits<Point_2>::Kernel Traits;
  typedef typename Traits::Segment_2 Segment_2;

public:
  PolylinesGraphicsItem(P* p_);

  void modelChanged();

public:
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  

  const QPen& verticesPen() const
  {
    return vertices_pen;
  }

  void setVerticesPen(const QPen& pen)
  {
    vertices_pen = pen;
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

  const QPen& edgesPen() const
  {
    return edges_pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    edges_pen = pen;
  }

protected:
  void updateBoundingBox();

  P * polylines;
  QPainter* m_painter;
  PainterOstream<Traits> painterostream;


  QRectF bounding_rect;

  QPen vertices_pen;
  QPen edges_pen;
  bool draw_edges;
  bool draw_vertices;
};


template <typename P>
PolylinesGraphicsItem<P>::PolylinesGraphicsItem(P * p_)
  :  polylines(p_), painterostream(0),
     draw_edges(true), draw_vertices(true)   
{
  setVerticesPen(QPen(::Qt::red, 3.));
  updateBoundingBox();
  setZValue(3);
}

template <typename P>
QRectF 
PolylinesGraphicsItem<P>::boundingRect() const
{
  return bounding_rect;
}




template <typename P>
void 
PolylinesGraphicsItem<P>::paint(QPainter *painter, 
                                const QStyleOptionGraphicsItem * /*option*/,
                                QWidget * /*widget*/)
{

  painterostream = PainterOstream<Traits>(painter);

  painter->setPen(this->edgesPen());  
  for(typename P::iterator it = polylines->begin();
      it != polylines->end();
      it++){
    Polyline & pl = *it;
    typename Polyline::iterator pit = pl.begin();
    Point_2 p = *pit;
    ++pit;
    for(; pit != pl.end(); ++pit){
      const Point_2 &q = *pit;
      painterostream << Segment_2(p, q);
      p = q;
    }
  }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename P>
void 
PolylinesGraphicsItem<P>::updateBoundingBox()
{
  Converter<Traits> convert;
  prepareGeometryChange();
  if(polylines->size() == 0){
    return;
  }
  Bbox_2 bb = CGAL::bounding_box(polylines->begin()->begin(), polylines->begin()->end()).bbox();
  for(typename P::iterator it = polylines->begin();
      it != polylines->end();
      ++it){
    bb = bb + CGAL::bounding_box(it->begin(), it->end()).bbox();;
  }

  bounding_rect = convert(bb);
}


template <typename P>
void 
PolylinesGraphicsItem<P>::modelChanged()
{
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_POLYLINES_GRAPHICS_ITEM_H
