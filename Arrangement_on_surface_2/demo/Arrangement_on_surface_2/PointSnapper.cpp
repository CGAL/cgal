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

#include <CGAL/ipower.h>
#include <QGraphicsView>

#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"
#include "PointSnapper.h"

#include <boost/optional.hpp>

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

PointSnapperBase::PointSnapperBase(
  QGraphicsScene* scene, GridGraphicsItem* grid) :
    GraphicsSceneMixin(scene),
    gridGraphicsItem{grid}, snapToGridEnabled{false},
    snapToArrangementEnabled{false}
{
}

namespace
{
struct ExplicitLambda
{
  template <typename Arrangement>
  void operator()(demo_types::TypeHolder<Arrangement>)
  {
    Arrangement* arr = nullptr;
    CGAL::assign(arr, arr_obj);
    res = new PointSnapper<Arrangement>(scene, grid, arr);
  }

  PointSnapperBase*& res;
  QGraphicsScene* scene;
  GridGraphicsItem* grid;
  CGAL::Object& arr_obj;
};
} // anonymous namespace

PointSnapperBase* PointSnapperBase::create(
  demo_types::TraitsType tt, QGraphicsScene* scene, GridGraphicsItem* grid,
  CGAL::Object arr_obj)
{
  PointSnapperBase* res;
  ExplicitLambda explicit_lambda{res, scene, grid, arr_obj};
  demo_types::visitArrangementType(tt, explicit_lambda);
  return res;
}

auto PointSnapperBase::snapPoint(const QPointF& qpt) -> Point_2
{
  if (this->snapToGridEnabled && this->snapToArrangementEnabled)
  {
    Point_2 pt = {qpt.x(), qpt.y()};
    auto gridPt = snapToGrid(qpt);
    auto arrPt = snapToArrangement(qpt);
    if (!arrPt)
      return gridPt;
    else if (
      compute_squared_distance_2(pt, gridPt) <
      compute_squared_distance_2(pt, *arrPt))
      return gridPt;
    else
      return *arrPt;
  }
  else if (this->snapToGridEnabled)
  {
    return snapToGrid(qpt);
  }
  else if (this->snapToArrangementEnabled)
  {
    auto arrPt = snapToArrangement(qpt);
    if (arrPt)
      return *snapToArrangement(qpt);
    else
      return {qpt.x(), qpt.y()};
  }
  else
    return Point_2{qpt.x(), qpt.y()};
}

// snap to grid without loss of precision
auto PointSnapperBase::snapToGrid(const QPointF& qpt) -> Point_2
{
  Rational two{2};
  Rational five{5};
  Rational half{0.5};
  // can't use 0.2 since it's not perfectly representable as a float/double
  Rational fifth{Rational{1} / five};

  Rational x;
  {
    int a = gridGraphicsItem->getXPower2();
    int b = gridGraphicsItem->getXPower5();
    // we have to calculate l in BigRat to be exact
    Rational lx;
    if (a < 0)
      lx = CGAL::ipower(half, -a);
    else
      lx = CGAL::ipower(two, a);
    if (b < 0)
      lx *= CGAL::ipower(fifth, -b);
    else
      lx *= CGAL::ipower(five, b);
    x = lx * std::lround(CGAL::to_double(Rational(qpt.x()) / lx));
  }

  Rational y;
  {
    int a = gridGraphicsItem->getYPower2();
    int b = gridGraphicsItem->getYPower5();
    // we have to calculate l in BigRat to be exact
    Rational ly;
    if (a < 0)
      ly = CGAL::ipower(half, -a);
    else
      ly = CGAL::ipower(two, a);
    if (b < 0)
      ly *= CGAL::ipower(fifth, -b);
    else
      ly *= CGAL::ipower(five, b);
    y = ly * std::lround(CGAL::to_double(Rational(qpt.y()) / ly));
  }

  return Point_2{x, y};
}

void PointSnapperBase::setSnapToGrid(bool val)
{
  this->snapToGridEnabled = val;
}
void PointSnapperBase::setSnapToArrangement(bool val)
{
  this->snapToArrangementEnabled = val;
}
bool PointSnapperBase::isSnapToGridEnabled() { return this->snapToGridEnabled; }
bool PointSnapperBase::isSnapToArrangementEnabled()
{
  return this->snapToArrangementEnabled;
}

template <typename Arr_>
PointSnapper<Arr_>::PointSnapper(
  QGraphicsScene* scene, GridGraphicsItem* grid, Arrangement* arr_) :
    PointSnapperBase(scene, grid),
    arr{arr_}
{
}

template <typename Arrangement>
inline boost::optional<PointSnapperBase::Point_2> snapToArrangement(
  const QPointF& qpt, const QTransform& worldTransform, Arrangement* arr)
{
  using Point_2 = PointSnapperBase::Point_2;
  using Rational = PointSnapperBase::Rational;
  using Compute_squared_distance_2 =
    PointSnapperBase::Compute_squared_distance_2;

  Compute_squared_distance_2 compute_squared_distance_2;

  Point_2 initialPoint{qpt.x(), qpt.y()};
  Point_2 closestPoint = initialPoint;

  bool first = true;
  Rational minDist(0);

  static constexpr int PIXEL_RADIUS = 15;
  Rational maxDist =
    PIXEL_RADIUS * PIXEL_RADIUS *
    std::abs(1 / (worldTransform.m11() * worldTransform.m22()));

  for (auto vit = arr->vertices_begin(); vit != arr->vertices_end(); ++vit)
  {
    auto arr_point = vit->point();
    Point_2 point{arr_point.x(), arr_point.y()};
    auto dist = compute_squared_distance_2(initialPoint, point);
    if (first || (dist < minDist))
    {
      first = false;
      minDist = dist;
      closestPoint = point;
    }
  }

  if (!first && minDist < maxDist)
    return closestPoint;
  else
    return {};
}

template <typename Traits>
struct SnapToArrangement
{
  using Point_2 = PointSnapperBase::Point_2;
  template <typename Arrangement>
  boost::optional<Point_2>
  operator()(const QPointF& qpt, const QTransform&, Arrangement*)
  {
    return Point_2{qpt.x(), qpt.y()};
  }
};
template <typename Kernel>
struct SnapToArrangement<CGAL::Arr_linear_traits_2<Kernel>>
{
  using Point_2 = PointSnapperBase::Point_2;
  template <typename Arrangement>
  boost::optional<Point_2> operator()(
    const QPointF& qpt, const QTransform& worldTransform, Arrangement* arr)
  {
    return snapToArrangement(qpt, worldTransform, arr);
  }
};
template <typename Kernel>
struct SnapToArrangement<CGAL::Arr_segment_traits_2<Kernel>>
{
  using Point_2 = PointSnapperBase::Point_2;
  template <typename Arrangement>
  boost::optional<Point_2> operator()(
    const QPointF& qpt, const QTransform& worldTransform, Arrangement* arr)
  {
    return snapToArrangement(qpt, worldTransform, arr);
  }
};
template <typename Kernel>
struct SnapToArrangement<CGAL::Arr_polyline_traits_2<Kernel>>
{
  using Point_2 = PointSnapperBase::Point_2;
  template <typename Arrangement>
  boost::optional<Point_2> operator()(
    const QPointF& qpt, const QTransform& worldTransform, Arrangement* arr)
  {
    return snapToArrangement(qpt, worldTransform, arr);
  }
};

template <typename Arr_>
auto PointSnapper<Arr_>::snapToArrangement(const QPointF& qpt)
  -> boost::optional<Point_2>
{
  using Traits = typename Arrangement::Geometry_traits_2;
  auto view = getView();
  if (view)
    return SnapToArrangement<Traits>{}(qpt, view->transform(), arr);
  else
    return {};
}
