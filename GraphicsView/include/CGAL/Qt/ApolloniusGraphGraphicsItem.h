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

#ifndef CGAL_QT_APOLLONIUS_GRAPH_GRAPHICS_ITEM_H
#define CGAL_QT_APOLLONIUS_GRAPH_GRAPHICS_ITEM_H

#include <CGAL/license/GraphicsView.h>




#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/utility.h>

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>

class QGraphicsSceneMouseEvent;


namespace CGAL {
namespace Qt {

template <typename AG, typename K>
class ApolloniusGraphGraphicsItem : public GraphicsItem
{
public:
  ApolloniusGraphGraphicsItem(AG* ag);


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

  const QPen& sitesPen() const
  {
    return sites_pen;
  }

  void setSitesPen(const QPen& pen)
  {
    sites_pen = pen;
  }

private:
  AG * ag;
  QPen edges_pen, sites_pen;
};



template <typename AG, typename K>
ApolloniusGraphGraphicsItem<AG,K>::ApolloniusGraphGraphicsItem(AG * ag)
  :  ag(ag)
{
  setZValue(3);
}

template <typename AG, typename K>
QRectF 
ApolloniusGraphGraphicsItem<AG,K>::boundingRect() const
{
  QRectF rect = CGAL::Qt::viewportsBbox(scene());
  return rect;
}


template <typename AG, typename K>
void 
ApolloniusGraphGraphicsItem<AG,K>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget * /*w*/)
{

  painter->setPen(this->sitesPen());
  QRectF rect = option->exposedRect;
  PainterOstream<K> pos(painter, rect);
  for(typename AG::Sites_iterator it = ag->sites_begin();
	 it != ag->sites_end(); it++ ) {
    pos << typename K::Circle_2(it->point(), square( it->weight()));
  }

  painter->setPen(this->edgesPen());
  ag->draw_dual(pos);
}


  template <typename AG, typename K>
void 
  ApolloniusGraphGraphicsItem<AG,K>::modelChanged()
{
  update();
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_APOLLONIUS_GRAPH_GRAPHICS_ITEM_H
