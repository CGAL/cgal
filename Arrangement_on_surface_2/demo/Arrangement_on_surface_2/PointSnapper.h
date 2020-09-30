// Copyright (c) 2020 GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Ahmed Essam <theartful.ae@gmail.com>

#ifndef ARRANGEMENT_DEMO_POINT_SNAPPER_H
#define ARRANGEMENT_DEMO_POINT_SNAPPER_H

#include "ArrangementTypes.h"
#include "GraphicsSceneMixin.h"
#include "GridGraphicsItem.h"

#include <boost/optional.hpp>

class GridGraphicsItem;
class QGraphicsScene;

class PointSnapperBase : public GraphicsSceneMixin
{
public:
  using Rational = demo_types::Rational;
  using Rat_kernel = demo_types::Rat_kernel;
  using Point_2 = demo_types::Rat_point_2;
  using Compute_squared_distance_2 = Rat_kernel::Compute_squared_distance_2;

public:
  PointSnapperBase(QGraphicsScene* scene, GridGraphicsItem* grid);

  Point_2 snapPoint(const QPointF& qpt);
  Point_2 snapToGrid(const QPointF& qpt);
  virtual boost::optional<Point_2> snapToArrangement(const QPointF& qpt) = 0;
  void setSnapToGrid(bool val);
  void setSnapToArrangement(bool val);
  bool isSnapToGridEnabled();
  bool isSnapToArrangementEnabled();

protected:
  GridGraphicsItem* gridGraphicsItem;
  bool snapToGridEnabled;
  bool snapToArrangementEnabled;
  Compute_squared_distance_2 compute_squared_distance_2;
};

template <typename Arr_>
class PointSnapper : public PointSnapperBase
{
  using Arrangement = Arr_;
  using Compute_squared_distance_2 = Rat_kernel::Compute_squared_distance_2;

public:
  PointSnapper(QGraphicsScene*, GridGraphicsItem*, Arrangement*);
  boost::optional<Point_2> snapToArrangement(const QPointF& qpt) override;

private:
  Arrangement* arr;
};

#endif
