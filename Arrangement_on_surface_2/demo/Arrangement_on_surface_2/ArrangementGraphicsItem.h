// Copyright (c) 2008, 2012, 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
#define CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H

#include <CGAL/Object.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <QImage>
#include <QPen>

#include "GraphicsSceneMixin.h"


namespace demo_types
{
enum class TraitsType : int;
}

namespace CGAL
{
namespace Qt
{

class ArrangementGraphicsItemBase :
    public GraphicsItem,
    public GraphicsSceneMixin
{
public:
  static ArrangementGraphicsItemBase*
  create(demo_types::TraitsType, CGAL::Object arr);

  const QPen& getVerticesPen() const;
  const QPen& getEdgesPen() const;
  void setVerticesPen(const QPen& pen);
  void setEdgesPen(const QPen& pen);
  void setBackgroundColor(QColor color);

  virtual QRectF getInterestingViewport() const = 0;

protected:
  ArrangementGraphicsItemBase();

  QPen verticesPen;
  QPen edgesPen;
  QPen facesPen;
}; // class ArrangementGraphicsItemBase

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_GRAPHICS_ITEM_H
