// Copyright (c) 2008  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ophir Setter <ophirset@post.tau.ac.il>
//

#ifndef CGAL_QT_SET_GRAPHICS_ITEM_H
#define CGAL_QT_SET_GRAPHICS_ITEM_H



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

template <typename Set>
class SetGraphicsItem : public GraphicsItem
{
public:
  SetGraphicsItem(const Set *arr);

  QRectF
  boundingRect() const;

  void
  paint(QPainter *painter,
        const QStyleOptionGraphicsItem *option,
        QWidget *widget);

  void
  modelChanged();

  const QPen& pen() const {
    return m_pen;
  }

  void setPen(const QPen& pen) {
    m_pen = pen;
  }

private:
  const Set *m_set;
  QPen m_pen;
};

template <typename Set>
SetGraphicsItem<Set>::SetGraphicsItem(const Set *set)
  : m_set(set)
{
  setZValue(3);
}

template <typename Set>
QRectF
SetGraphicsItem<Set>::boundingRect() const
{
  QRectF rect = CGAL::Qt::viewportsBbox(scene());
  return rect;
}


template <typename Set>
void
SetGraphicsItem<Set>::paint(QPainter *painter,
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget* )
{
  if (m_set == nullptr)
    return;

  QRectF rect = option->exposedRect;
  // R is the kernel. Move it to a template parameter.
  PainterOstream<typename std::iterator_traits<typename Set::const_iterator>::value_type::R>
    pos(painter, rect);

  painter->setPen(this->pen());
  for(typename Set::const_iterator it = m_set->begin();
      it != m_set->end(); it++) {
    pos << *it;
  }
}


template <typename T>
void
SetGraphicsItem<T>::modelChanged()
{
  update();
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_SET_GRAPHICS_ITEM_H
