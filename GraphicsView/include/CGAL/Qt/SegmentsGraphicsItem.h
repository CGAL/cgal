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

#ifndef CGAL_QT_SEGMENTS_GRAPHICS_ITEM_H
#define CGAL_QT_SEGMENTS_GRAPHICS_ITEM_H

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
class SegmentsGraphicsItem : public GraphicsItem
{
  typedef typename P::value_type Segment_2;
  typedef typename CGAL::Kernel_traits<Segment_2>::Kernel Traits;

public:
  SegmentsGraphicsItem(P* p_);

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

protected:
  void updateBoundingBox();

  P * segments;
  QPainter* m_painter;
  PainterOstream<Traits> painterostream;


  QRectF bounding_rect;

  QPen vertices_pen;
  bool draw_edges;
  bool draw_vertices;
};


template <typename P>
SegmentsGraphicsItem<P>::SegmentsGraphicsItem(P * p_)
  :  segments(p_), painterostream(0),
     draw_edges(true), draw_vertices(true)   
{
  setVerticesPen(QPen(::Qt::red, 3.));

  updateBoundingBox();
  setZValue(3);
}

template <typename P>
QRectF 
SegmentsGraphicsItem<P>::boundingRect() const
{
  return bounding_rect;
}




template <typename P>
void 
SegmentsGraphicsItem<P>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem * /*option*/,
                                    QWidget * /*widget*/)
{

  painterostream = PainterOstream<Traits>(painter);
  painter->setPen(QPen(::Qt::black, 0, ::Qt::SolidLine, ::Qt::RoundCap, ::Qt::RoundJoin));
    for(typename P::iterator it = segments->begin();
        it != segments->end();
        it++){
      painterostream << *it;
    }
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename P>
void 
SegmentsGraphicsItem<P>::updateBoundingBox()
{
  Converter<Traits> convert;
  prepareGeometryChange();
  if(segments->size() == 0){
    return;
  }
  Bbox_2 bb = segments->begin()->bbox();
  for(typename P::iterator it = segments->begin();
      it != segments->end();
      ++it){
    bb = bb + it->bbox();
  }

  bounding_rect = convert(bb);
}


template <typename P>
void 
SegmentsGraphicsItem<P>::modelChanged()
{
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_SEGMENTS_GRAPHICS_ITEM_H
