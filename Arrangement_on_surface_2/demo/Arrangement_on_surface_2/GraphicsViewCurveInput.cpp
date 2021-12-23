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
#include "AlgebraicCurveParser.h"
#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"
#include "CurveInputMethods.h"
#include "GraphicsViewCurveInputTyped.h"
#include "QtMetaTypes.h"
#include "Utils/Utils.h"

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

template <typename Arr_>
GraphicsViewCurveInput<Arr_>::GraphicsViewCurveInput(
  Arrangement* arrangement_, QObject* parent, QGraphicsScene* scene) :
    GraphicsViewCurveInputBase(parent, scene),
    arrangement(arrangement_)
{
  this->setDefaultInputMethod(
    std::integral_constant<
      bool, std::tuple_size<InputMethodTuple>::value != 0>{});
  for_each(inputMethods, [&](auto&& it) {
    it.setScene(scene);
    it.setCallback(this);
  });
  curveGenerator.setTraits(this->arrangement->traits());
}

template <typename Arr_>
void GraphicsViewCurveInput<Arr_>::setCurveType(CurveType type)
{
  this->reset();
  for_each(inputMethods, [&](auto&& it) {
    if (it.curveType() == type)
      this->setInputMethod(static_cast<CurveInputMethod*>(&it));
  });
}

template <typename Arr_>
void GraphicsViewCurveInput<Arr_>::setPointSnapper(PointSnapperBase* snapper_)
{
  for_each(inputMethods, [&](auto&& it) { it.setPointSnapper(snapper_); });
}

template <typename Arr_>
template <typename>
void GraphicsViewCurveInput<Arr_>::setDefaultInputMethod(std::true_type)
{
  this->setInputMethod(&std::get<0>(inputMethods));
}

template <typename Arr_>
void GraphicsViewCurveInput<Arr_>::setDefaultInputMethod(std::false_type)
{
}

template <typename Arr_>
void GraphicsViewCurveInput<Arr_>::generate(CGAL::Object o)
{
  insertCurve(
    demo_types::enumFromArrType<Arrangement>(),
    CGAL::make_object(this->arrangement), o);
  Q_EMIT CGAL::Qt::GraphicsViewCurveInputBase::modelChanged();
}

template <typename Arr_>
void GraphicsViewCurveInput<Arr_>::curveInputDoneEvent(
  const std::vector<Point_2>& clickedPoints, CurveType type)
{
  boost::optional<Curve_2> cv =
    this->curveGenerator.generate(clickedPoints, type);
  if (cv)
  {
    Insert_curve<Arrangement>{}(this->arrangement, *cv);
    Q_EMIT this->modelChanged();
  }
}

// CurveGeneratorBase
template <typename ArrTraits_>
auto CurveGeneratorBase<ArrTraits_>::generate(
  const std::vector<Point_2>& clickedPoints, CurveType type)
  -> boost::optional<Curve_2>
{
  boost::optional<Curve_2> res;
  switch (type)
  {
  case CurveType::Segment:
    res = generateSegment(clickedPoints);
    break;
  case CurveType::Ray:
    res = generateRay(clickedPoints);
    break;
  case CurveType::Line:
    res = generateLine(clickedPoints);
    break;
  case CurveType::Polyline:
    res = generatePolyline(clickedPoints);
    break;
  case CurveType::Circle:
    res = generateCircle(clickedPoints);
    break;
  case CurveType::Ellipse:
    res = generateEllipse(clickedPoints);
    break;
  case CurveType::ThreePointCircularArc:
    res = generateThreePointCircularArc(clickedPoints);
    break;
  case CurveType::FivePointConicArc:
    res = generateFivePointConicArc(clickedPoints);
    break;
  case CurveType::Bezier:
    res = generateBezier(clickedPoints);
  default:
    break;
  }
  return res;
}

template <typename ArrTraits_>
void CurveGeneratorBase<ArrTraits_>::setTraits(const ArrTraits* traits_)
{
  this->traits = traits_;
}

// Curve Generator Segment Traits
template <typename Kernel_>
auto CurveGenerator<CGAL::Arr_segment_traits_2<Kernel_>>::generateSegment(
  const std::vector<Point_2>& clickedPoints) -> boost::optional<Curve_2>
{
  Curve_2 res{clickedPoints[0], clickedPoints[1]};
  return res;
}

// Curve Generator Polyline Traits
template <typename SegmentTraits>
auto CurveGenerator<CGAL::Arr_polyline_traits_2<SegmentTraits>>::
  generatePolyline(const std::vector<Point_2>& clickedPoints)
    -> boost::optional<Curve_2>
{
  if (clickedPoints.size() < 2) return {};

  auto construct_poly = this->traits->construct_curve_2_object();
  Curve_2 res = construct_poly(clickedPoints.begin(), clickedPoints.end());
  return res;
}

// Curve Generator Linear Traits
template <typename Kernel_>
auto CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateSegment(
  const std::vector<Point_2>& points) -> boost::optional<Curve_2>
{
  Curve_2 res = Curve_2(Segment_2(points[0], points[1]));
  return res;
}

template <typename Kernel_>
auto CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateRay(
  const std::vector<Point_2>& points) -> boost::optional<Curve_2>
{
  Curve_2 res = Curve_2(Ray_2(points[0], points[1]));
  return res;
}

template <typename Kernel_>
auto CurveGenerator<CGAL::Arr_linear_traits_2<Kernel_>>::generateLine(
  const std::vector<Point_2>& points) -> boost::optional<Curve_2>
{
  Curve_2 res = Curve_2(Line_2(points[0], points[1]));
  return res;
}

// CurveGenerator Conic Traits
template <typename RatKernel, typename AlgKernel, typename NtTraits>
auto CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateSegment(const std::vector<Point_2>& points)
    -> boost::optional<Curve_2>
{
  Curve_2 res = Curve_2(Rat_segment_2(points[0], points[1]));
  return res;
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
auto CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateCircle(const std::vector<Point_2>& points) -> boost::optional<Curve_2>
{
  auto sq_rad =
    (points[0].x() - points[1].x()) * (points[0].x() - points[1].x()) +
    (points[0].y() - points[1].y()) * (points[0].y() - points[1].y());
  Curve_2 res = Curve_2(Rat_circle_2(points[0], sq_rad));
  return res;
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
auto CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateEllipse(const std::vector<Point_2>& points)
    -> boost::optional<Curve_2>
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
  return res;
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
auto CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateThreePointCircularArc(const std::vector<Point_2>& points)
    -> boost::optional<Curve_2>
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
    return res;
  }
  else
  {
    std::cout << "Points don't specify a valid conic." << std::endl;
    return {};
  }
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
auto CurveGenerator<Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
  generateFivePointConicArc(const std::vector<Point_2>& points)
    -> boost::optional<Curve_2>
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
      return res;
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
auto CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateLine(const std::vector<Point_2>& points) -> boost::optional<Curve_2>
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
  return res;
}

template <typename Coefficient_>
auto CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateCircle(const std::vector<Point_2>& points) -> boost::optional<Curve_2>
{
  auto sq_rad =
    (points[0].x() - points[1].x()) * (points[0].x() - points[1].x()) +
    (points[0].y() - points[1].y()) * (points[0].y() - points[1].y());
  return this->generateEllipse_(points[0], sq_rad, sq_rad);
}

template <typename Coefficient_>
auto CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateEllipse(const std::vector<Point_2>& points)
    -> boost::optional<Curve_2>
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
auto CurveGenerator<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  generateEllipse_(const Point_2& center, Rational rxRat, Rational ryRat)
    -> boost::optional<Curve_2>
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
  return res;
}

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
auto CurveGenerator<
  Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits, BoundingTraits>>::
  generateBezier(const std::vector<Point_2>& clickedPoints)
    -> boost::optional<Curve_2>
{
  if (clickedPoints.size() < 2) return {};
  return Curve_2{clickedPoints.begin(), clickedPoints.end()};
}

// msvc2015 doesn't play well with polymorphic lambdas
namespace
{
struct ExplicitLambda
{
  template <typename Arrangement>
  void operator()(demo_types::TypeHolder<Arrangement>)
  {
    Arrangement* arr = nullptr;
    CGAL::assign(arr, arr_obj);
    res = new GraphicsViewCurveInput<Arrangement>(arr, parent, scene);
  }

  GraphicsViewCurveInputBase*& res;
  CGAL::Object& arr_obj;
  QObject* parent;
  QGraphicsScene* scene;
};
} // anonymous namespace

GraphicsViewCurveInputBase* GraphicsViewCurveInputBase::create(
  demo_types::TraitsType tt, CGAL::Object arr_obj, QObject* parent,
  QGraphicsScene* scene)
{
  GraphicsViewCurveInputBase* res;
  ExplicitLambda explicit_lambda{res, arr_obj, parent, scene};
  demo_types::visitArrangementType(tt, explicit_lambda);
  return res;
}

} // namespace Qt
} // namespace CGAL

#ifdef CGAL_USE_CORE
CGAL::Object algebraicCurveFromExpression(
  const CGAL::Object& arr_obj, const std::string& exp, bool& is_first_curve)
{
  using Polynomial_2 = demo_types::DemoTypes::Alg_seg_traits::Polynomial_2;
  using Alg_seg_arr = demo_types::DemoTypes::Alg_seg_arr;

  Alg_seg_arr* arr;
  if (!CGAL::assign(arr, arr_obj)) CGAL_error();

  is_first_curve = (arr->number_of_edges() == 0);

  auto poly = AlgebraicCurveParser<Polynomial_2>{}(exp);
  if (!poly) return {};

  auto construct_curve = arr->traits()->construct_curve_2_object();
  auto cv = construct_curve(*poly);

  return CGAL::make_object(cv);
}

CGAL::Object rationalCurveFromExpression(
  const CGAL::Object& arr_obj, const std::string& numerator,
  const std::string& denominator, bool& is_first_curve)
{
  using Polynomial_1 = demo_types::DemoTypes::Rational_traits::Polynomial_1;
  using Rational_arr = demo_types::DemoTypes::Rational_arr;

  Rational_arr* arr;
  if (!CGAL::assign(arr, arr_obj)) CGAL_error();

  is_first_curve = (arr->number_of_edges() == 0);

  auto curve_parser = AlgebraicCurveParser<Polynomial_1>{};
  auto poly_num = curve_parser(numerator);
  if (!poly_num) return {};
  auto poly_den = curve_parser(denominator);
  if (!poly_den) return {};

  auto construct_curve = arr->traits()->construct_curve_2_object();
  auto cv = construct_curve(*poly_num, *poly_den);

  return CGAL::make_object(cv);
}
#endif
