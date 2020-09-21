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

#include "SplitEdgeCallback.h"
#include "ArrangementTypes.h"
#include "PointSnapper.h"
#include "Utils.h"

#include <CGAL/polynomial_utils.h>

#include <QGraphicsSceneMouseEvent>
#include <QGraphicsScene>


SplitEdgeCallbackBase::SplitEdgeCallbackBase( QObject* parent ) :
  CGAL::Qt::Callback( parent )
{
  this->setColor(::Qt::darkGray);
  segmentInputMethod.setCallback(this);
}

bool SplitEdgeCallbackBase::eventFilter(QObject* object, QEvent* event)
{
  return segmentInputMethod.eventFilter(object, event);
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
  intersectCurves( this->arr->traits()->intersect_2_object( ) ),
  areEqual( this->arr->traits()->equal_2_object( ) )
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
void SplitEdgeCallback<Arr_>::curveInputDoneEvent(
  const std::vector<Input_point_2>& clickedPoints, CGAL::Qt::CurveType)
{
  try
  {
    this->splitEdges(clickedPoints[0], clickedPoints[1]);
  }
  catch (const std::exception& ex)
  {
    std::cerr << ex.what() << '\n';
    std::cerr << __FILE__ << ':' << __LINE__ << '\n';
    return;
  }
}

template <typename Traits_>
struct ConstructSegment
{
  using Traits = Traits_;

  template <typename Point_2>
  auto operator()(const Traits* traits, const Point_2& p1, const Point_2& p2)
  {
    Arr_construct_point_2<Traits> toArrPoint{traits};

    auto construct_x_monotone_curve_2 =
      traits->construct_x_monotone_curve_2_object();
    return construct_x_monotone_curve_2(toArrPoint(p1), toArrPoint(p2));
  }
};

template <typename Coefficient_>
struct ConstructSegment<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>
{
  using Traits = CGAL::Arr_algebraic_segment_traits_2<Coefficient_>;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;

  template <typename Point_2>
  auto operator()(const Traits* traits, const Point_2& p1, const Point_2& p2)
  {
    Arr_construct_point_2<Traits> toArrPoint{traits};

    auto construct_x_monotone_segment_2 =
      traits->construct_x_monotone_segment_2_object();

    std::vector<X_monotone_curve_2> curves;
    construct_x_monotone_segment_2(
      toArrPoint(p1), toArrPoint(p2), std::back_inserter(curves));
    return curves[0];
  }
};

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
struct ConstructSegment<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>>
{
  using Traits = CGAL::Arr_Bezier_curve_traits_2<
    RatKernel, AlgKernel, NtTraits, BoundingTraits>;
  using Curve_2 = typename Traits::Curve_2;
  using X_monotone_curve_2 = typename Traits::X_monotone_curve_2;

  template <typename Point_2>
  auto operator()(const Traits* traits, const Point_2& p1, const Point_2& p2)
  {
    Arr_construct_point_2<Traits> toArrPoint{traits};

    typename Traits::Point_2 points[] = {toArrPoint(p1), toArrPoint(p2)};
    Curve_2 curve{points, points + 2};
    auto make_x_monotone = traits->make_x_monotone_2_object();
    std::vector<CGAL::Object> curves;
    make_x_monotone(curve, std::back_inserter(curves));
    X_monotone_curve_2 segment;
    CGAL::assign(segment, curves[0]);
    return segment;
  }
};

template <typename AlgebraicKernel_d_1>
struct ConstructSegment<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>
{
  using Traits = CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>;
  using Rational = typename Traits::Rational;
  using Integer = typename Traits::Integer;
  using Polynomial_1 = typename Traits::Polynomial_1;
  using Algebraic_real_1 = typename Traits::Algebraic_real_1;
  using RationalTraits = CGAL::Rational_traits<typename Traits::Rational>;

  template <typename Point_2>
  auto operator()(const Traits* traits, const Point_2& p1, const Point_2& p2)
  {
    RationalTraits ratTraits;

    Rational dx = p2.x() - p1.x();
    Rational dy = p2.y() - p1.y();
    Polynomial_1 x = CGAL::shift(Polynomial_1(1), 1, 0);
    Polynomial_1 poly_num;
    Polynomial_1 poly_den;

    if (dx != 0)
    {
      Rational mRat = dy / dx;
      Rational cRat = p1.y() - mRat * p1.x();
      // y = (a/b) x + (e/f)
      auto a = ratTraits.numerator(mRat);
      auto b = ratTraits.denominator(mRat);
      auto e = ratTraits.numerator(cRat);
      auto f = ratTraits.denominator(cRat);

      poly_num = f * a * x + b * e;
      poly_den = Polynomial_1{b * f};
    }
    else
    {
      throw std::runtime_error(
        "Vertical split segments are not allowed with rational traits!\n");
    }

    auto construct_x_monotone_curve_2 =
      traits->construct_x_monotone_curve_2_object();

    return construct_x_monotone_curve_2(
      poly_num, poly_den, Algebraic_real_1{p1.x()}, Algebraic_real_1{p2.x()});
  }
};

template <typename Arr_>
void SplitEdgeCallback<Arr_>::splitEdges(
  const Input_point_2& p1, const Input_point_2& p2)
{
  X_monotone_curve_2 splitCurve =
    ConstructSegment<Traits>{}(this->arr->traits(), p1, p2);

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

ARRANGEMENT_DEMO_SPECIALIZE_ARR(SplitEdgeCallback)
