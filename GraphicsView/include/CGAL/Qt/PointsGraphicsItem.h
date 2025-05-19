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

#ifndef CGAL_QT_POINTS_GRAPHICS_ITEM_H
#define CGAL_QT_POINTS_GRAPHICS_ITEM_H

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
class PointsGraphicsItem : public GraphicsItem
{
  typedef typename std::iterator_traits<typename P::iterator>::value_type Point_2;
  typedef typename CGAL::Kernel_traits<Point_2>::Kernel Traits;

public:
  PointsGraphicsItem(P* p_);

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

  P * points;
  QPainter* m_painter;
  PainterOstream<Traits> painterostream;


  QRectF bounding_rect;

  QPen vertices_pen;
  bool draw_vertices;
};


template <typename P>
PointsGraphicsItem<P>::PointsGraphicsItem(P * p_)
  :  points(p_), painterostream(0),  draw_vertices(true)
{
  setVerticesPen(QPen(::Qt::red, 3.));
  if(points->size() == 0){
    this->hide();
  }
  updateBoundingBox();
  setZValue(3);
}

template <typename P>
QRectF
PointsGraphicsItem<P>::boundingRect() const
{
  return bounding_rect;
}




template <typename P>
void
PointsGraphicsItem<P>::paint(QPainter *painter,
                                    const QStyleOptionGraphicsItem * /*option*/,
                                    QWidget * /*widget*/)
{
  if(drawVertices()) {
    Converter<Traits> convert;

    painter->setPen(verticesPen());
    QTransform matrix = painter->worldTransform();
    painter->resetTransform();
    for(typename P::iterator it = points->begin();
        it != points->end();
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
PointsGraphicsItem<P>::updateBoundingBox()
{
  Converter<Traits> convert;
  prepareGeometryChange();
  if(points->size() == 0){
    return;
  }
  bounding_rect = convert(CGAL::bounding_box(points->begin(), points->end()));
}


template <typename P>
void
PointsGraphicsItem<P>::modelChanged()
{
  if((points->size() == 0) ){
    this->hide();
  } else if((points->size() > 0) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_POINTS_GRAPHICS_ITEM_H
