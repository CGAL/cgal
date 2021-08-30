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

#include "CurveGraphicsItem.h"
#include "ArrangementPainterOstream.h"
#include "ArrangementTypes.h"
#include "Utils/ConstructBoundingBox.h"
#include "Utils/Utils.h"

#include <limits>

namespace CGAL
{
namespace Qt
{

template <class ArrTraits>
CurveGraphicsItem<ArrTraits>::CurveGraphicsItem() :
    bb(), m_edgeColor(::Qt::red), m_edgeWidth(2),
    m_vertexColor(::Qt::red), m_vertexRadius(1)
{
  this->setZValue(4);
  this->pointsGraphicsItem.setParentItem(this);
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::paint(
  QPainter* painter, const QStyleOptionGraphicsItem* /* option */,
  QWidget* /* widget */)
{
  // draw the curves
  QPen edgesPen(this->m_edgeColor, this->m_edgeWidth);
  edgesPen.setCosmetic(true);
  painter->setPen(edgesPen);

  auto painterOstream = ArrangementPainterOstream<Traits>(painter);
  painterOstream.setScene(this->getScene());

  for (auto& curve : this->curves) { painterOstream << curve; }
}

template <class ArrTraits>
QRectF CurveGraphicsItem<ArrTraits>::boundingRect() const
{
  auto viewport = this->viewportRect();
  qreal xmin = viewport.left();
  qreal ymin = viewport.top();
  qreal xmax = viewport.right();
  qreal ymax = viewport.bottom();
  if (this->bb.xmin() > xmin) xmin = this->bb.xmin();
  if (this->bb.ymin() > ymin) ymin = this->bb.ymin();
  if (this->bb.xmax() < xmax) xmax = this->bb.xmax();
  if (this->bb.ymax() < ymax) ymax = this->bb.ymax();
  if (xmin > xmax || ymin > ymax)
  {
    xmin = 0;
    xmax = 0;
    ymin = 0;
    ymax = 0;
  }
  return {QPointF{xmin, ymin}, QPointF{xmax, ymax}};
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::insert(const X_monotone_curve_2& curve)
{
  this->curves.push_back(curve);
  this->updateBoundingBox();
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::insert(const Point_2& point)
{
  this->pointsGraphicsItem.insert(point);
  this->updateBoundingBox();
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::clear()
{
  this->curves.clear();
  this->pointsGraphicsItem.clear();

  this->updateBoundingBox();
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::modelChanged()
{
  this->updateBoundingBox();
  this->update();
}

template <class ArrTraits>
const QColor& CurveGraphicsItem<ArrTraits>::edgeColor() const
{
  return this->m_edgeColor;
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::setEdgeColor(const QColor& color)
{
  this->m_edgeColor = color;
}

template <class ArrTraits>
int CurveGraphicsItem<ArrTraits>::edgeWidth() const
{
  return this->m_edgeWidth;
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::setEdgeWidth(int width)
{
  this->m_edgeWidth = width;
}

template <class ArrTraits>
const QColor& CurveGraphicsItem<ArrTraits>::vertexColor() const
{
  return this->m_vertexColor;
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::setVertexColor(const QColor& color)
{
  this->m_vertexColor = color;
}

template <class ArrTraits>
int CurveGraphicsItem<ArrTraits>::vertexRadius() const
{
  return this->m_vertexRadius;
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::setVertexRadius(int radius)
{
  this->m_vertexRadius = radius;
}

template <class ArrTraits>
void CurveGraphicsItem<ArrTraits>::updateBoundingBox()
{
  this->prepareGeometryChange();

  this->bb = {};
  ConstructBoundingBox<Traits> construct_bounding_box;
  for (auto& curve : curves)
    this->bb += construct_bounding_box(curve);
}

ARRANGEMENT_DEMO_SPECIALIZE_TRAITS(CurveGraphicsItem)

} // namespace Qt
} // namespace CGAL
