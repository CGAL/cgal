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

namespace CGAL {
namespace Qt {

//
template <typename GeometryTraits>
CurveGraphicsItem<GeometryTraits>::CurveGraphicsItem(const Traits& traits) :
  m_traits(traits),
  bb(),
  m_edgeColor(::Qt::red),
  m_edgeWidth(2),
  m_vertexColor(::Qt::red),
  m_vertexRadius(1)
{
  this->setZValue(4);
  this->pointsGraphicsItem.setParentItem(this);
}

//
template <typename GeometryTraits>
void CurveGraphicsItem<GeometryTraits>::
paint(QPainter* painter,
      const QStyleOptionGraphicsItem* /* option */,
      QWidget* /* widget */) {
  // draw the curves
  QPen edgesPen(this->m_edgeColor, this->m_edgeWidth);
  edgesPen.setCosmetic(true);
  painter->setPen(edgesPen);

  auto painterOstream = ArrangementPainterOstream<Traits>(painter);
  painterOstream.setScene(this->getScene());

  for (auto& curve : this->curves) { painterOstream << curve; }
}

//
template <typename GeometryTraits>
QRectF CurveGraphicsItem<GeometryTraits>::boundingRect() const {
  auto viewport = this->viewportRect();
  qreal xmin = viewport.left();
  qreal ymin = viewport.top();
  qreal xmax = viewport.right();
  qreal ymax = viewport.bottom();
  if (this->bb.xmin() > xmin) xmin = this->bb.xmin();
  if (this->bb.ymin() > ymin) ymin = this->bb.ymin();
  if (this->bb.xmax() < xmax) xmax = this->bb.xmax();
  if (this->bb.ymax() < ymax) ymax = this->bb.ymax();
  if (xmin > xmax || ymin > ymax) {
    xmin = 0;
    xmax = 0;
    ymin = 0;
    ymax = 0;
  }
  return {QPointF{xmin, ymin}, QPointF{xmax, ymax}};
}

//
template <typename GeometryTraits>
void CurveGraphicsItem<GeometryTraits>::insert(const X_monotone_curve_2& curve) {
  this->curves.push_back(curve);
  this->updateBoundingBox();
}

//
template <typename GeometryTraits>
void CurveGraphicsItem<GeometryTraits>::insert(const Point_2& point) {
  this->pointsGraphicsItem.insert(point);
  this->updateBoundingBox();
}

//
template <typename GeometryTraits>
void CurveGraphicsItem<GeometryTraits>::clear() {
  this->curves.clear();
  this->pointsGraphicsItem.clear();

  this->updateBoundingBox();
}

//
template <typename GeometryTraits>
void CurveGraphicsItem<GeometryTraits>::modelChanged() {
  this->updateBoundingBox();
  this->update();
}

//
template <typename GeometryTraits>
const QColor& CurveGraphicsItem<GeometryTraits>::edgeColor() const {
  return this->m_edgeColor;
}

//
template <typename GeometryTraits>
void CurveGraphicsItem<GeometryTraits>::setEdgeColor(const QColor& color) {
  this->m_edgeColor = color;
}

//
template <typename GeometryTraits>
int CurveGraphicsItem<GeometryTraits>::edgeWidth() const {
  return this->m_edgeWidth;
}

//
template <typename GeometryTraits>
void CurveGraphicsItem<GeometryTraits>::setEdgeWidth(int width) {
  this->m_edgeWidth = width;
}

//
template <typename GeometryTraits>
const QColor& CurveGraphicsItem<GeometryTraits>::vertexColor() const {
  return this->m_vertexColor;
}

//
template <typename GeometryTraits>
void CurveGraphicsItem<GeometryTraits>::setVertexColor(const QColor& color) {
  this->m_vertexColor = color;
}

//
template <typename GeometryTraits>
int CurveGraphicsItem<GeometryTraits>::vertexRadius() const {
  return this->m_vertexRadius;
}

//
template <typename GeometryTraits>
void CurveGraphicsItem<GeometryTraits>::setVertexRadius(int radius) {
  this->m_vertexRadius = radius;
}

//
template <typename GeometryTraits>
void CurveGraphicsItem<GeometryTraits>::updateBoundingBox() {
  this->prepareGeometryChange();

  this->bb = {};
  ConstructBoundingBox<Traits> ctr_bbox(m_traits);
  for (auto& curve : curves) this->bb += ctr_bbox(curve);
}

ARRANGEMENT_DEMO_SPECIALIZE_TRAITS(CurveGraphicsItem)

} // namespace Qt
} // namespace CGAL
