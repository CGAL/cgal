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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_VORONOI_GRAPHICS_ITEM_H
#define CGAL_QT_VORONOI_GRAPHICS_ITEM_H



#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/utility.h>

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>

#include <CGAL/intersection_2.h>

class QGraphicsSceneMouseEvent;


namespace CGAL {
namespace Qt {

template <typename DT>
class VoronoiGraphicsItem : public GraphicsItem
{
public:
  VoronoiGraphicsItem(DT  * dt_);


  QRectF 
  boundingRect() const;
  
  void 
  paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  
  void 
  modelChanged();

  const QPen& edgesPen() const
  {
    return edges_pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    edges_pen = pen;
  }

private:

  void updateCCenters() {
    cc.clear();
    for (typename DT::Face_iterator it = dt->faces_begin(); it != dt->faces_end(); it++) {
      cc[it] = Circumcenter()(*it);
    }
  }

  typedef typename DT::Hyperbolic_Voronoi_point                                   Voronoi_point;
  typedef typename DT::Geom_traits::Construct_inexact_hyperbolic_circumcenter_2   Circumcenter;
  typedef typename DT::Geom_traits::Construct_hyperbolic_segment_2                Segment;
  typedef typename DT::Face_handle                                                Face_handle;
  typedef typename DT::Geom_traits::Construct_hyperbolic_point_2                  CP2;

  DT * dt;
  QPen edges_pen;
  std::map<Face_handle, Voronoi_point> cc;   // to hold the circumcenters of all the faces
};



template <typename DT>
VoronoiGraphicsItem<DT>::VoronoiGraphicsItem(DT * dt_)
  :  dt(dt_)
{
  setZValue(3);
  updateCCenters();
}

template <typename DT>
QRectF 
VoronoiGraphicsItem<DT>::boundingRect() const
{
  QRectF rect = CGAL::Qt::viewportsBbox(scene());
  return rect;
}


template <typename DT>
void 
VoronoiGraphicsItem<DT>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget * /*w*/)
{
  QRectF rect = option->exposedRect;
  PainterOstream<typename DT::Geom_traits> pos(painter, rect);
  
  painter->setPen(edgesPen());
  
  // delete
  QPen temp = painter->pen();
  QPen old = temp;
  temp.setWidthF(0.01);
  painter->setPen(temp);
  //
  
  for (typename DT::Face_iterator fit = dt->faces_begin(); fit != dt->faces_end(); fit++) {
    for (int i = 0; i < 3; i++) {
      typename DT::Hyperbolic_segment s =  Segment()(cc[fit], CP2()(cc[fit->neighbor(i)], dt->neighbor_translation(fit, i))); //dt->dual(std::pair<typename DT::Face_handle, int>(fit, i));
      pos << s;
    }
  }

  painter->setPen(old);
}


template <typename T>
void 
VoronoiGraphicsItem<T>::modelChanged()
{
  updateCCenters();
  update();
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_VORONOI_GRAPHICS_ITEM_H
