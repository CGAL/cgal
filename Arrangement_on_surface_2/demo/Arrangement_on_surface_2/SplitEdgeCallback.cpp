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

#include "SplitEdgeCallback.h"
#include "Utils.h"
#include "PointSnapper.h"

#include <QGraphicsSceneMouseEvent>
#include <QGraphicsScene>


SplitEdgeCallbackBase::SplitEdgeCallbackBase( QObject* parent ) :
  CGAL::Qt::Callback( parent )
{
  this->setColor(::Qt::darkGray);
  segmentInputMethod.setCallback(this);
}

void SplitEdgeCallbackBase::setColor(QColor c)
{
  this->segmentInputMethod.setColor(c);
}

QColor SplitEdgeCallbackBase::getColor( ) const
{
  return this->segmentInputMethod.getColor();
}


void SplitEdgeCallbackBase::setPointSnapper(PointSnapperBase* snapper_)
{
  this->segmentInputMethod.setPointSnapper(snapper_);
}

template <typename Arr_>
SplitEdgeCallback<Arr_>::SplitEdgeCallback(Arrangement* arr_, QObject* parent):
  SplitEdgeCallbackBase( parent ),
  arr( arr_ ),
  traits( arr->geometry_traits() ),
  intersectCurves( this->traits->intersect_2_object( ) ),
  areEqual( this->traits->equal_2_object( ) )
{
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::setScene( QGraphicsScene* scene_ )
{
  CGAL::Qt::Callback::setScene(scene_);
  this->segmentInputMethod.setScene(scene_);
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::reset( )
{
  this->segmentInputMethod.reset();
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
  this->segmentInputMethod.mouseMoveEvent(event);
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
  this->segmentInputMethod.mousePressEvent(event);
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::curveInputDoneEvent(
  const std::vector<Input_point_2>& clickedPoints, CGAL::Qt::CurveType)
{
  typedef typename ArrTraitsAdaptor<Traits>::CoordinateType CoordinateType;

  auto pt1 = clickedPoints[0];
  auto pt2 = clickedPoints[1];

  this->splitEdges(
    {CoordinateType{pt1.x()}, CoordinateType{pt1.y()}},
    {CoordinateType{pt2.x()}, CoordinateType{pt2.y()}}, traits);
}

template <typename Traits_>
struct ConstructSegment
{
  using Traits = Traits_;
  using Point_2 = typename Traits::Point_2;

  auto operator()(const Traits* traits, const Point_2& p1, const Point_2& p2)
  {
    auto construct_x_monotone_curve_2 =
      traits->construct_x_monotone_curve_2_object();
    return construct_x_monotone_curve_2(p1, p2);
  }
};

template <typename RatKernel, class AlgKernel, class NtTraits>
struct ConstructSegment<
  CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>>
{
  using Traits =
    CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>;
  using Curve_2 = typename Traits::Curve_2;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;
  using Point_2 = typename Traits::Point_2;

  auto operator()(const Traits* traits, const Point_2& p1, const Point_2& p2)
  {
    Point_2 points[] = {p1, p2};
    Curve_2 curve{points, points + 2};
    auto make_x_monotone = traits->make_x_monotone_2_object();
    std::vector<CGAL::Object> curves;
    make_x_monotone(curve, std::back_inserter(curves));
    X_monotone_curve_2 segment;
    CGAL::assign(segment, curves[0]);
    return segment;
  }
};

template <typename Arr_>
template <typename TTraits>
void SplitEdgeCallback<Arr_>::splitEdges(
  const Point_2& p1, const Point_2& p2, const TTraits*)
{
  X_monotone_curve_2 splitCurve = ConstructSegment<TTraits>{}(traits, p1, p2);

  for (auto hei = this->arr->halfedges_begin();
       hei != this->arr->halfedges_end(); ++hei)
  {
    X_monotone_curve_2 curve = hei->curve();
    CGAL::Object res;
    CGAL::Oneset_iterator<CGAL::Object> oi(res);
    this->intersectCurves(splitCurve, curve, oi);
    std::pair<Point_2, Multiplicity> pair;
    if (hei == this->arr->halfedges_end()) continue;
    if (CGAL::assign(pair, res))
    {
      Point_2 splitPoint = pair.first;
      if (
        (!hei->source()->is_at_open_boundary() &&
         this->areEqual(hei->source()->point(), splitPoint)) ||
        (!hei->target()->is_at_open_boundary() &&
         this->areEqual(hei->target()->point(), splitPoint)))
      { continue; }
      this->arr->split_edge(hei, splitPoint);
    }
  }
  this->reset();
  Q_EMIT modelChanged();
}

template <typename Arr_>
template <typename Coefficient_>
void SplitEdgeCallback<Arr_>::
splitEdges(const Point_2& p1, const Point_2& p2,
           const CGAL::Arr_algebraic_segment_traits_2<Coefficient_>*)
{
  typename Traits::Construct_x_monotone_segment_2 constructSegment =
    traits->construct_x_monotone_segment_2_object();

  std::vector<X_monotone_curve_2> curves;

  constructSegment(p1, p2, std::back_inserter(curves));

  X_monotone_curve_2 splitCurve = curves[0];

  for (auto hei = this->arr->halfedges_begin();
       hei != this->arr->halfedges_end(); ++hei)
  {
    X_monotone_curve_2 curve = hei->curve();
    CGAL::Object res;
    CGAL::Oneset_iterator<CGAL::Object> oi(res);
    this->intersectCurves(splitCurve, curve, oi);
    std::pair<Point_2, Multiplicity> pair;
    if (CGAL::assign(pair, res))
    {
      Point_2 splitPoint = pair.first;
      if (
        (!hei->source()->is_at_open_boundary() &&
         this->areEqual(hei->source()->point(), splitPoint)) ||
        (!hei->target()->is_at_open_boundary() &&
         this->areEqual(hei->target()->point(), splitPoint)))
      { continue; }
      this->arr->split_edge(hei, splitPoint);
    }
  }
  this->reset();
  Q_EMIT modelChanged();
}

ARRANGEMENT_DEMO_SPECIALIZE_ARR(SplitEdgeCallback)
