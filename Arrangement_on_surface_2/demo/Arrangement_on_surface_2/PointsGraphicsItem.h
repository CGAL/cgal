// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef POINTS_GRAPHICS_ITEM_H
#define POINTS_GRAPHICS_ITEM_H

#include <vector>
#include <QGraphicsItem>
#include <CGAL/number_utils.h>

class QPainter;

/**
   Add a set of points to the QGraphicsScene.
*/
class PointsGraphicsItem : public QGraphicsItem
{
public:
  PointsGraphicsItem(QGraphicsItem* parent = nullptr);

  void paint(
    QPainter* painter, const QStyleOptionGraphicsItem* option = nullptr,
    QWidget* widget = nullptr) override;

  QRectF boundingRect() const override;

  template <class Point>
  inline void insert(const Point& point);

  void clear();

  void setColor(QColor c);
  QColor getColor() const;

  void setPointRadius(double d);
  double getPointRadius() const;

protected:
  std::vector<QPointF> points;
  double pointRadius;
  QColor color;

}; // class PointsGraphicsItem

namespace CGAL
{
template <
  class RatKernel_, class AlgKernel_, class NtTraits_, class BoundingTraits_>
class _Bezier_point_2;
}

namespace details
{
struct ApproximatePoint
{
  template <typename T>
  QPointF operator()(const T& point)
  {
    return QPointF(CGAL::to_double(point.x()), CGAL::to_double(point.y()));
  }

  template <
    class RatKernel, class AlgKernel, class NtTraits, class BoundingTraits>
  QPointF operator()(
    const CGAL::_Bezier_point_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>&
      point)
  {
    auto xy = point.approximate();
    return QPointF(xy.first, xy.second);
  }
};
}

template <class Point>
inline void PointsGraphicsItem::insert(const Point& point)
{
  this->prepareGeometryChange();
  this->points.push_back(details::ApproximatePoint{}(point));
}

#endif // POINTS_GRAPHICS_ITEM_H
