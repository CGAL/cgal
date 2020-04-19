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

#ifndef CGAL_QT_POWERDIAGRAM_GRAPHICS_ITEM_H
#define CGAL_QT_POWERDIAGRAM_GRAPHICS_ITEM_H

#include <CGAL/license/GraphicsView.h>




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

template <typename RT>
class PowerdiagramGraphicsItem : public GraphicsItem
{

  typedef typename RT::Geom_traits Geom_traits;
  typedef typename Kernel_traits<typename RT::Bare_point>::Kernel K;
  typedef typename K::Segment_2 Segment_2;
  typedef typename K::Line_2 Line_2;
  typedef typename K::Ray_2 Ray_2;

public:
  PowerdiagramGraphicsItem(RT  * rt_);


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
  RT * rt;
  QPen edges_pen;
};



template <typename RT>
PowerdiagramGraphicsItem<RT>::PowerdiagramGraphicsItem(RT * rt_)
  :  rt(rt_)
{
  setZValue(3);
}

template <typename RT>
QRectF
PowerdiagramGraphicsItem<RT>::boundingRect() const
{
  QRectF rect = CGAL::Qt::viewportsBbox(scene());

  return rect;
}


template <typename RT>
void
PowerdiagramGraphicsItem<RT>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget * /*w*/)
{
  QRectF rect = option->exposedRect;
  PainterOstream<K> pos(painter, rect);

  painter->setPen(edgesPen());
  for(typename RT::Finite_edges_iterator eit = rt->finite_edges_begin();
      eit != rt->finite_edges_end();
      eit++){
    CGAL::Object o = rt->dual(eit);
    Segment_2 s;
    Ray_2 r;
    Line_2 l;
    if(CGAL::assign(s,o)){
      pos << s;
    } else if(CGAL::assign(r,o)) {
      pos << r;
    }else if(CGAL::assign(l,o)) {
      pos << l;
    }
  }
}


template <typename T>
void
PowerdiagramGraphicsItem<T>::modelChanged()
{
  update();
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_POWERDIAGRAM_GRAPHICS_ITEM_H
