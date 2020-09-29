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

#include "GraphicsViewCurveInput.h"
#include "ArrangementTypes.h"
#include "QtMetaTypes.h"

#include <CGAL/polynomial_utils.h>

#include <QEvent>
#include <QKeyEvent>

template <std::size_t I = 0, typename FuncT, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
for_each(std::tuple<Tp...>&, FuncT)
{
}

template <std::size_t I = 0, typename FuncT, typename... Tp>
  inline typename std::enable_if <
  I<sizeof...(Tp), void>::type for_each(std::tuple<Tp...>& t, FuncT f)
{
  f(std::get<I>(t));
  for_each<I + 1, FuncT, Tp...>(t, f);
}

namespace CGAL
{
namespace Qt
{

GraphicsViewCurveInputBase::GraphicsViewCurveInputBase(
  QObject* parent, QGraphicsScene* scene) :
    Callback(parent, scene),
    inputMethod(nullptr)
{
}

void GraphicsViewCurveInputBase::setInputMethod(CurveInputMethod* inputMethod_)
{
  this->inputMethod = inputMethod_;
}

void GraphicsViewCurveInputBase::reset()
{
  if (this->inputMethod)
  {
    this->inputMethod->reset();
    this->inputMethod = nullptr;
  }
}

bool GraphicsViewCurveInputBase::eventFilter(QObject* obj, QEvent* event)
{
  return this->inputMethod->eventFilter(obj, event);
}

void GraphicsViewCurveInputBase::setColor(QColor c)
{
  this->inputMethod->setColor(c);
}

template <typename ArrTraits>
GraphicsViewCurveInput<ArrTraits>::GraphicsViewCurveInput(
  const ArrTraits* traits, QObject* parent, QGraphicsScene* scene) :
    GraphicsViewCurveInputBase(parent, scene)
{
  this->setDefaultInputMethod(
    std::integral_constant<
      bool, std::tuple_size<InputMethodTuple>::value != 0>{});
  QObject::connect(
    &curveGenerator, SIGNAL(generate(CGAL::Object)), this,
    SIGNAL(generate(CGAL::Object)));
  for_each(inputMethods, [&](auto&& it) {
    it.setScene(scene);
    it.setCallback(&curveGenerator);
  });
  curveGenerator.setTraits(traits);
}

template <typename ArrTraits>
void GraphicsViewCurveInput<ArrTraits>::setCurveType(CurveType type)
{
  this->reset();
  for_each(inputMethods, [&](auto&& it) {
    if (it.curveType() == type)
      this->setInputMethod(static_cast<CurveInputMethod*>(&it));
  });
}

template <typename ArrTraits>
void GraphicsViewCurveInput<ArrTraits>::setPointSnapper(
  PointSnapperBase* snapper_)
{
  for_each(inputMethods, [&](auto&& it) { it.setPointSnapper(snapper_); });
}

template <typename ArrTraits>
template <typename>
void GraphicsViewCurveInput<ArrTraits>::setDefaultInputMethod(std::true_type)
{
  this->setInputMethod(&std::get<0>(inputMethods));
}

template <typename ArrTraits>
void GraphicsViewCurveInput<ArrTraits>::setDefaultInputMethod(std::false_type)
{
}

// CurveGeneratorBase
void CurveGeneratorBase::curveInputDoneEvent(
  const std::vector<Point_2>& clickedPoints, CurveType type)
{
  CGAL::Object obj;
  switch (type)
  {
  case CurveType::Segment:
    obj = generateSegment(clickedPoints);
    break;
  case CurveType::Ray:
    obj = generateRay(clickedPoints);
    break;
  case CurveType::Line:
    obj = generateLine(clickedPoints);
    break;
  case CurveType::Polyline:
    obj = generatePolyline(clickedPoints);
    break;
  case CurveType::Circle:
    obj = generateCircle(clickedPoints);
    break;
  case CurveType::Ellipse:
    obj = generateEllipse(clickedPoints);
    break;
  case CurveType::ThreePointCircularArc:
    obj = generateThreePointCircularArc(clickedPoints);
    break;
  case CurveType::FivePointConicArc:
    obj = generateFivePointConicArc(clickedPoints);
    break;
  case CurveType::Bezier:
    obj = generateBezier(clickedPoints);
  default:
    break;
  }
  if (!obj.empty()) Q_EMIT generate(obj);
}

CGAL::Object CurveGeneratorBase::generateSegment(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generatePolyline(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generateRay(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generateLine(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generateCircle(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generateEllipse(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object
CurveGeneratorBase::generateThreePointCircularArc(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object
CurveGeneratorBase::generateFivePointConicArc(const std::vector<Point_2>&)
{
  return {};
}
CGAL::Object CurveGeneratorBase::generateBezier(const std::vector<Point_2>&)
{
  return {};
}

template <typename ArrTraits_>
void CurveGeneratorTypedBase<ArrTraits_>::setTraits(const ArrTraits* traits_)
{
  this->traits = traits_;
}

// Curve Generator Segment Traits
template <typename Kernel_>
CGAL::Object
CurveGenerator<CGAL::Arr_segment_traits_2<Kernel_>>::generateSegment(
  const std::vector<Point_2>& clickedPoints)
{
  Curve_2 res{clickedPoints[0], clickedPoints[1]};
  return CGAL::make_object(res);
}

// Curve Generator Polyline Traits
template <typename SegmentTraits>
CGAL::Object
CurveGenerator<CGAL::Arr_polyline_traits_2<SegmentTraits>>::generatePolyline(
  const std::vector<Point_2>& clickedPoints)
{
  if (clickedPoints.size() < 2) return {};

  auto construct_poly = this->traits->construct_curve_2_object();
  Curve_2 res = construct_poly(clickedPoints.begin(), clickedPoints.end());
  return CGAL::make_object(res);
}

// Curve Generator Linear Traits
template <typename Kernel_>
CGAL::Object
CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateSegment(
  const std::vector<Point_2>& points)
{
  Curve_2 res = Curve_2(Segment_2(points[0], points[1]));
  return CGAL::make_object(res);
}

template <typename Kernel_>
CGAL::Object CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateRay(
  const std::vector<Point_2>& points)
{
  Curve_2 res = Curve_2(Ray_2(points[0], points[1]));
  return CGAL::make_object(res);
}

template <typename Kernel_>
CGAL::Object CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateLine(
  const std::vector<Point_2>& points)
{
  Curve_2 res = Curve_2(Line_2(points[0], points[1]));
  return CGAL::make_object(res);
}

// CurveGenerator Conic Traits
template <typename RatKernel, typename AlgKernel, typename NtTraits>
CGAL::Object
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateSegment(const std::vector<Point_2>& points)
{
  Curve_2 res = Curve_2(Rat_segment_2(points[0], points[1]));
  return CGAL::make_object(res);
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
CGAL::Object
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateCircle(const std::vector<Point_2>& points)
{
  auto sq_rad =
    (points[0].x() - points[1].x()) * (points[0].x() - points[1].x()) +
    (points[0].y() - points[1].y()) * (points[0].y() - points[1].y());
  Curve_2 res = Curve_2(Rat_circle_2(points[0], sq_rad));
  return CGAL::make_object(res);
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
CGAL::Object
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateEllipse(const std::vector<Point_2>& points)
{
  auto x1 = (CGAL::min)(points[0].x(), points[1].x());
  auto y1 = (CGAL::min)(points[0].y(), points[1].y());
  auto x2 = (CGAL::max)(points[0].x(), points[1].x());
  auto y2 = (CGAL::max)(points[0].y(), points[1].y());

  Rat_FT a = CGAL::abs(Rat_FT(x1) - Rat_FT(x2)) / 2;
  Rat_FT b = CGAL::abs(Rat_FT(y1) - Rat_FT(y2)) / 2;
  Rat_FT a_sq = a * a;
  Rat_FT b_sq = b * b;
  Rat_FT x0 = (x2 + x1) / 2;
  Rat_FT y0 = (y2 + y1) / 2;

  Rat_FT r = b_sq;
  Rat_FT s = a_sq;
  Rat_FT t = 0;
  Rat_FT u = -2 * x0 * b_sq;
  Rat_FT v = -2 * y0 * a_sq;
  Rat_FT ww = x0 * x0 * b_sq + y0 * y0 * a_sq - a_sq * b_sq;

  Curve_2 res = Curve_2(r, s, t, u, v, ww);
  return CGAL::make_object(res);
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
CGAL::Object
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateThreePointCircularArc(const std::vector<Point_2>& points)
{
  auto& qp1 = points[0];
  auto& qp2 = points[1];
  auto& qp3 = points[2];
  Rat_point_2 p1 = Rat_point_2(qp1.x(), qp1.y());
  Rat_point_2 p2 = Rat_point_2(qp2.x(), qp2.y());
  Rat_point_2 p3 = Rat_point_2(qp3.x(), qp3.y());
  RatKernel ker;
  if (!ker.collinear_2_object()(p1, p2, p3))
  {
    Curve_2 res(p1, p2, p3);
    return CGAL::make_object(res);
  }
  else
  {
    std::cout << "Points don't specify a valid conic." << std::endl;
    return {};
  }
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
CGAL::Object
CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateFivePointConicArc(const std::vector<Point_2>& points)
{
  auto& qp0 = points[0];
  auto& qp1 = points[1];
  auto& qp2 = points[2];
  auto& qp3 = points[3];
  auto& qp4 = points[4];
  Rat_point_2 p0 = Rat_point_2(qp0.x(), qp0.y());
  Rat_point_2 p1 = Rat_point_2(qp1.x(), qp1.y());
  Rat_point_2 p2 = Rat_point_2(qp2.x(), qp2.y());
  Rat_point_2 p3 = Rat_point_2(qp3.x(), qp3.y());
  Rat_point_2 p4 = Rat_point_2(qp4.x(), qp4.y());
  try
  {
    Curve_2 res(p0, p1, p2, p3, p4);
    if (res.is_valid())
      return CGAL::make_object(res);
    else
      std::cout << "Points don't specify a valid conic. Try again!"
                << std::endl;
  }
  catch (...)
  {
    std::cout << "Points don't specify a valid conic. Try again!" << std::endl;
  }
  return {};
}

// CurveGenerator Algebraic Traits
template <typename Coefficient_>
CGAL::Object
CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateLine(const std::vector<Point_2>& points)
{
  RationalTraits ratTraits;

  Rational dx = points[1].x() - points[0].x();
  Rational dy = points[1].y() - points[0].y();
  Polynomial_2 x = CGAL::shift(Polynomial_2(1), 1, 0);
  Polynomial_2 y = CGAL::shift(Polynomial_2(1), 1, 1);

  Polynomial_2 poly;
  if (dx != 0)
  {
    Rational mRat = dy / dx;
    Rational cRat = points[0].y() - mRat * points[0].x();
    // y = (a/b) x + (e/f)
    auto a = ratTraits.numerator(mRat);
    auto b = ratTraits.denominator(mRat);
    auto e = ratTraits.numerator(cRat);
    auto f = ratTraits.denominator(cRat);

    poly = b * f * y - f * a * x - b * e;
  }
  // vertical line
  else
  {
    Rational xP = points[0].x();
    auto a = ratTraits.numerator(xP);
    auto b = ratTraits.denominator(xP);
    poly = b * x - a;
  }

  auto construct_curve = this->traits->construct_curve_2_object();
  auto res = construct_curve(poly);
  return CGAL::make_object(res);
}

template <typename Coefficient_>
CGAL::Object
CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateCircle(const std::vector<Point_2>& points)
{
  auto sq_rad =
    (points[0].x() - points[1].x()) * (points[0].x() - points[1].x()) +
    (points[0].y() - points[1].y()) * (points[0].y() - points[1].y());
  return this->generateEllipse_(points[0], sq_rad, sq_rad);
}

template <typename Coefficient_>
CGAL::Object
CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateEllipse(const std::vector<Point_2>& points)
{
  auto rx =
    (points[0].x() - points[1].x()) * (points[0].x() - points[1].x()) / 4.;
  auto ry =
    (points[0].y() - points[1].y()) * (points[0].y() - points[1].y()) / 4.;
  Point_2 center = {
    (points[1].x() + points[0].x()) / 2, (points[1].y() + points[0].y()) / 2};
  return this->generateEllipse_(center, rx, ry);
}

template <typename Coefficient_>
CGAL::Object
CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateEllipse_(const Point_2& center, Rational rxRat, Rational ryRat)
{
  RationalTraits ratTraits;

  Polynomial_2 x = CGAL::shift(Polynomial_2(1), 1, 0);
  Polynomial_2 y = CGAL::shift(Polynomial_2(1), 1, 1);
  Rational xCenter = center.x();
  Rational yCenter = center.y();

  // (g/h) (x - a/b)^2 + (e/f) (y - c/d)^2 = (e/f) * (g/h)
  auto a = ratTraits.numerator(xCenter);
  auto b = ratTraits.denominator(xCenter);
  auto c = ratTraits.numerator(yCenter);
  auto d = ratTraits.denominator(yCenter);
  auto e = ratTraits.numerator(rxRat);
  auto f = ratTraits.denominator(rxRat);
  auto g = ratTraits.numerator(ryRat);
  auto h = ratTraits.denominator(ryRat);

  Polynomial_2 poly = g * f * CGAL::ipower(b * d * x - d * a, 2) +
                      e * h * CGAL::ipower(b * d * y - b * c, 2) -
                      e * g * b * b * d * d;

  auto construct_curve = this->traits->construct_curve_2_object();
  auto res = construct_curve(poly);
  return CGAL::make_object(res);
}

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
CGAL::Object CurveGenerator<
  Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>>::
  generateBezier(const std::vector<Point_2>& clickedPoints)
{
  using Traits =
    Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>;

  if (clickedPoints.size() < 2) return {};

  return CGAL::make_object(
    typename Traits::Curve_2{clickedPoints.begin(), clickedPoints.end()});
}

ARRANGEMENT_DEMO_SPECIALIZE_TRAITS(GraphicsViewCurveInput)

} // namespace Qt
} // namespace CGAL
