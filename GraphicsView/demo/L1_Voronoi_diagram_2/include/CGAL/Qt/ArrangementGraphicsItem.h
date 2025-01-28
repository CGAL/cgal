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

#ifndef CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
#define CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H



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

template <typename Arr>
class ArrangementGraphicsItem : public GraphicsItem
{
public:
  ArrangementGraphicsItem(const Arr *arr);

  void setArrangement(const Arr *arr);

  QRectF
  boundingRect() const;

  void
  paint(QPainter *painter,
        const QStyleOptionGraphicsItem *option,
        QWidget *widget);

  void
  modelChanged();

  const QPen& edgesPen() const
  {
    return m_edges_pen;
  }

  const QPen& verticesPen() const
  {
    return m_vertices_pen;
  }

  void setEdgesPen(const QPen& pen)
  {
    m_edges_pen = pen;
  }

  void setVerticesPen(const QPen& pen)
  {
    m_vertices_pen = pen;
  }

private:
  const Arr *m_arr;
  QPen m_edges_pen;
  QPen m_vertices_pen;
};

template <typename Arr>
ArrangementGraphicsItem<Arr>::ArrangementGraphicsItem(const Arr *arr)
  : m_arr(arr)
{
  setZValue(3);
}

template <typename Arr>
void ArrangementGraphicsItem<Arr>::setArrangement(const Arr *arr) {
  m_arr = arr;
}


template <typename Arr>
QRectF
ArrangementGraphicsItem<Arr>::boundingRect() const
{
  QRectF rect = CGAL::Qt::viewportsBbox(scene());
  return rect;
}

template <typename Arr>
void
ArrangementGraphicsItem<Arr>::paint(QPainter *painter,
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget* )
{
  if (m_arr == nullptr)
    return;

  QRectF rect = option->exposedRect;
  PainterOstream<typename Arr::Geometry_traits_2> pos(painter, rect);

  painter->setPen(edgesPen());
  for(typename Arr::Edge_const_iterator eit = m_arr->edges_begin();
      eit != m_arr->edges_end(); eit++) {
    pos << eit->curve();
  }

  painter->setPen(verticesPen());
  for(typename Arr::Vertex_const_iterator vit = m_arr->vertices_begin();
      vit != m_arr->vertices_end(); vit++) {
    pos << vit->point();
  }
}


template <typename T>
void
ArrangementGraphicsItem<T>::modelChanged()
{
  update();
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
