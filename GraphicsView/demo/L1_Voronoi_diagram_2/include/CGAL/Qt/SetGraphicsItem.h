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
// $URL: svn+ssh://ophirset@scm.gforge.inria.fr/svn/cgal/trunk/GraphicsView/include/CGAL/Qt/ArrangementGraphicsItem.h $
// $Id: ArrangementGraphicsItem.h 45924 2008-10-01 09:19:53Z afabri $
// 
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
                                    QWidget *w)
{
  if (m_set == NULL)
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
