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

#ifndef CGAL_QT_CIRCULAR_ARC_GRAPHICS_ITEM_H
#define CGAL_QT_CIRCULAR_ARC_GRAPHICS_ITEM_H

#include <CGAL/license/GraphicsView.h>


#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename CK>
class CircularArcGraphicsItem : public GraphicsItem
{
  typedef typename CK::Point_2 Point_2;
  typedef typename CK::Circle_2 Circle_2;
  typedef typename CK::Circular_arc_2 Circular_arc_2;

public:
  CircularArcGraphicsItem();

  void modelChanged();

public:
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  

  const QPen& verticesPen() const
  {
    return this->vertices_pen;
  }

  const QPen& edgesPen() const
  {
    return edges_pen;
  }

  void setVerticesPen(const QPen& pen)
  {
    this->vertices_pen = pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    edges_pen = pen;
  }
  
  void setArc(const Circular_arc_2& a);

  Circular_arc_2 arc() const
  {
    return arc_;
  }

protected:
  void updateBoundingBox();

  QPainter* m_painter;
  PainterOstream<CK> painterostream;

  QRectF bounding_rect;

  QPen edges_pen;
  QPen vertices_pen;

  Circular_arc_2 arc_;
};


template <typename CK>
void 
CircularArcGraphicsItem<CK>::setArc(const Circular_arc_2& a)
{
  arc_ = a;
  updateBoundingBox();
  update();
}

template <typename CK>
CircularArcGraphicsItem<CK>::CircularArcGraphicsItem()
  : painterostream(0)
{
  this->hide();
  setZValue(3);
}

template <typename CK>
QRectF 
CircularArcGraphicsItem<CK>::boundingRect() const
{
  return bounding_rect;
}




template <typename CK>
void 
CircularArcGraphicsItem<CK>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem * /*option*/,
                                    QWidget * /*widget*/)
{
  painter->setPen(this->edgesPen());
  painterostream = PainterOstream<CK>(painter);
  
  painterostream << arc_;
}

template <typename CK>
void 
CircularArcGraphicsItem<CK>::updateBoundingBox()
{
  Converter<CK> convert;
  prepareGeometryChange();
  //bounding_rect = convert(arc_.supporting_circle().bbox());
  bounding_rect = convert(arc_.bbox());
}


template <typename CK>
void 
CircularArcGraphicsItem<CK>::modelChanged()
{
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_CIRCULAR_ARC_GRAPHICS_ITEM_H
