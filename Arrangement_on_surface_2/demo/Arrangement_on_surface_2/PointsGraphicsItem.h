// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef POINTS_GRAPHICS_ITEM_H
#define POINTS_GRAPHICS_ITEM_H

#include <vector>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/number_utils.h>
#include <QPen>

class QPainter;
class QPen;

/**
   Add a set of points to the QGraphicsScene.
*/
class PointsGraphicsItem: public CGAL::Qt::GraphicsItem
{
public:
  PointsGraphicsItem( );

  void paint(
    QPainter* painter, const QStyleOptionGraphicsItem* option = nullptr,
    QWidget* widget = nullptr) override;

  QRectF boundingRect( ) const override;				//!< virtual function for the bounding box

  /** Template type
     *  adds the points to the vector
     */
  template < class Point >
  inline void insert( const Point& point );

  void clear( );

  void setColor( QColor c );				//!< sets the color of the curve.
  QColor getColor( ) const;					//!< returns the color of the curve

  void setPointRadius( double d );			//!< sets the user defined radius of the curve
  double getPointRadius( ) const;			//!< returns the radius of the curve

public Q_SLOTS:
  virtual void modelChanged( );

protected:
  std::vector< QPointF > points;  		/*!< vector of points of the curve */
  double pointRadius;					/*!< area of the curve draw */
  QColor color;                       	/*!< QColor object for the curve */

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
