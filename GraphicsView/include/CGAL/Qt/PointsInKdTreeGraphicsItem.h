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

#ifndef CGAL_QT_POINTS_IN_KD_TREE_GRAPHICS_ITEM_H
#define CGAL_QT_POINTS_IN_KD_TREE_GRAPHICS_ITEM_H

#include <CGAL/license/GraphicsView.h>


#include <CGAL/Bbox_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/iterator.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename KdTree>
class PointsInKdTreeGraphicsItem : public GraphicsItem
{
  typedef typename std::iterator_traits<typename KdTree::iterator>::value_type Point_2;
  typedef typename CGAL::Kernel_traits<Point_2>::Kernel Traits;
  typedef typename Traits::Iso_rectangle_2 Iso_rectangle_2;

  typedef CGAL::Fuzzy_iso_box<typename KdTree::Traits> Fuzzy_iso_box;


  // Instead of first collecting points into a container, and then draw them
  // we use an output iterator that draws them on the fly
  template <typename K>
  class Draw
    : public CGAL::cpp98::iterator<std::output_iterator_tag, void, void, void, void> {
    QPainter* painter;
    QTransform* matrix;
    Converter<K> convert;
  public:
    Draw(QPainter* painter, QTransform* matrix)
      : painter(painter), matrix(matrix)
    {}

    Draw& operator=(const Point_2& p)
    {
      QPointF point = matrix->map(convert(p));
      painter->drawPoint(point);
      return *this;
    }

    Draw& operator++()
    {
      return *this;
    }

    Draw& operator*()
    {
      return *this;
    }


    Draw& operator++(int)
    {
      return *this;
    }

  };


public:
  PointsInKdTreeGraphicsItem(KdTree* p_);

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

protected:
  void updateBoundingBox();

  KdTree * kdtree;
  QPainter* m_painter;
  PainterOstream<Traits> painterostream;
  Converter<Traits> convert;

  QRectF bounding_rect;

  QPen vertices_pen;
  bool draw_vertices;
};


template <typename KdTree>
PointsInKdTreeGraphicsItem<KdTree>::PointsInKdTreeGraphicsItem(KdTree * p_)
  :  kdtree(p_), painterostream(0),  draw_vertices(true)
{
  setVerticesPen(QPen(::Qt::red, 3.));
  if(kdtree->size() == 0){
    this->hide();
  }
  updateBoundingBox();
  setZValue(3);
  setFlag(QGraphicsItem::ItemUsesExtendedStyleOption, true);
}

template <typename KdTree>
QRectF
PointsInKdTreeGraphicsItem<KdTree>::boundingRect() const
{
  return bounding_rect;
}






template <typename KdTree>
void
PointsInKdTreeGraphicsItem<KdTree>::paint(QPainter *painter,
                                    const QStyleOptionGraphicsItem *option,
                                    QWidget * /*widget*/)
{
  Iso_rectangle_2 isor = convert(option->exposedRect);
  Fuzzy_iso_box range(isor.vertex(0), isor.vertex(2));
  painter->setPen(verticesPen());
  QTransform matrix = painter->worldTransform();
  painter->resetTransform();
  Draw<Traits> draw(painter, &matrix);
  kdtree->search(draw, range);
}

// We let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
template <typename KdTree>
void
PointsInKdTreeGraphicsItem<KdTree>::updateBoundingBox()
{
  prepareGeometryChange();
  if(kdtree->size() == 0){
    return;
  }
  bounding_rect = convert(CGAL::bounding_box(kdtree->begin(), kdtree->end()));
}


template <typename KdTree>
void
PointsInKdTreeGraphicsItem<KdTree>::modelChanged()
{
  if((kdtree->size() == 0) ){
    this->hide();
  } else if((kdtree->size() > 0) && (! this->isVisible())){
    this->show();
  }
  updateBoundingBox();
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_POINTS_IN_KD_TREE_GRAPHICS_ITEM_H
