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

#include "PointsGraphicsItem.h"

#include <limits>
#include <QPen>
#include <QPainter>

//! A constructor.
PointsGraphicsItem::PointsGraphicsItem(QGraphicsItem* parent) :
    QGraphicsItem(parent), pointRadius(3.0), color(::Qt::blue)
{
}

void PointsGraphicsItem::paint(
  QPainter* painter, const QStyleOptionGraphicsItem* /* option */,
  QWidget* /* widget */)
{
  painter->save();

  painter->setBrush(QBrush(this->color));
  QTransform matrix = painter->worldTransform();
  painter->resetTransform();

  for (auto& point : this->points)
    painter->drawEllipse(
      matrix.map(point), this->pointRadius, this->pointRadius);

  painter->restore();
}

QRectF PointsGraphicsItem::boundingRect() const
{
  if (this->points.size() == 0) { return QRectF(); }
  double xmin = (std::numeric_limits<double>::max)();
  double xmax = -(std::numeric_limits<double>::max)();
  double ymin = (std::numeric_limits<double>::max)();
  double ymax = -(std::numeric_limits<double>::max)();
  for (unsigned int i = 0; i < this->points.size(); ++i)
  {
    QPointF pt = this->points[i];
    double x = pt.x();
    double y = pt.y();
    xmin = (std::min)(xmin, x);
    xmax = (std::max)(xmax, x);
    ymin = (std::min)(ymin, y);
    ymax = (std::max)(ymax, y);
  }
  QRectF res(QPointF(xmin, ymin), QPointF(xmax, ymax));
  res.adjust(-5, -5, 5, 5); // pad the borders a bit
  return res;
}

void PointsGraphicsItem::clear( )
{
  this->prepareGeometryChange( );

  this->points.clear( );
}

void PointsGraphicsItem::setColor( QColor c )
{
  this->color = c;
}

QColor PointsGraphicsItem::getColor( ) const
{
  return this->color;
}

void PointsGraphicsItem::setPointRadius( double d )
{
  this->pointRadius = d;
}

double PointsGraphicsItem::getPointRadius( ) const
{
  return this->pointRadius;
}
