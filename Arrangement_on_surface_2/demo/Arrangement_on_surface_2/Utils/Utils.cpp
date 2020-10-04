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

#include "ArrangementPainterOstream.h"
#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"

#include "PointLocationFunctions.h"
#include "Utils.h"

#include <QGraphicsView>
#include <QScrollBar>

#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/Arr_overlay_2.h>

template <typename Kernel_>
double
Compute_squared_distance_2<CGAL::Arr_segment_traits_2<Kernel_>>::operator()(
  const Point_2& p, const X_monotone_curve_2& c) const
{
  Point_2 p1 = c.source();
  Point_2 p2 = c.target();
  Segment_2 seg(p1, p2);

  return CGAL::to_double(CGAL::squared_distance(p, seg));
}

template <typename Kernel_>
double
Compute_squared_distance_2<CGAL::Arr_linear_traits_2<Kernel_>>::operator()(
  const Point_2& p, const X_monotone_curve_2& c) const
{
  Segment_2 seg;
  Ray_2 ray;
  Line_2 line;
  FT res;
  if (c.is_segment())
  {
    seg = c.segment();
    res = CGAL::squared_distance(p, seg);
  }
  else if (c.is_ray())
  {
    ray = c.ray();
    res = CGAL::squared_distance(p, ray);
  }
  else // ( c.is_line( ) )
  {
    line = c.line();
    res = CGAL::squared_distance(p, line);
  }
  return CGAL::to_double(res);
}

template <typename Kernel_>
double
Compute_squared_distance_2<CGAL::Arr_polyline_traits_2<Kernel_>>::operator()(
  const Point_2& p, const X_monotone_curve_2& c) const
{
  Seg_const_it seg_it_s = c.subcurves_begin();

  bool first = true;
  FT min_dist = 0;

  while (seg_it_s != c.subcurves_end())
  {
    Segment_2 seg = *seg_it_s;
    FT dist = CGAL::squared_distance(p, seg);

    if (first || dist < min_dist)
    {
      first = false;
      min_dist = dist;
    }
    seg_it_s++;
  }

  return CGAL::to_double(min_dist);
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
double Compute_squared_distance_2<
  CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
operator()(const Point_2& p, const X_monotone_curve_2& c) const
{
  // Get the co-ordinates of the curve's source and target.
  // double sx = CGAL::to_double( c.source( ).x( ) );
  // double sy = CGAL::to_double( c.source( ).y( ) );
  // double tx = CGAL::to_double( c.target( ).x( ) );
  // double ty = CGAL::to_double( c.target( ).y( ) );

  if (c.orientation() == CGAL::COLLINEAR)
  {
    Point_2 ps = c.source();
    Point_2 pt = c.target();
    Segment_2 seg(ps, pt);

    FT res = CGAL::squared_distance(p, seg);
    return CGAL::to_double(res);
  }
  else
  {
    // If the curve is monotone, than its source and its target has the
    // extreme x co-ordinates on this curve.
    // bool is_source_left = (sx < tx);
    // int  x_min = is_source_left ? (*w).x_pixel(sx) : (*w).x_pixel(tx);
    // int  x_max = is_source_left ? (*w).x_pixel(tx) : (*w).x_pixel(sx);
    // double   prev_x = is_source_left ? sx : tx;
    // double   prev_y = is_source_left ? sy : ty;
    // double   curr_x, curr_y;
    // int      x;
    // Arr_conic_point_2 px;

    bool first = true;
    FT min_dist(100000000);
    // AlgKernel ker;

    int n = 100;
    if (this->scene != NULL && this->scene->views().size() != 0)
    { // use the scene to approximate the resolution of the curve
      QGraphicsView* view = this->scene->views().first();
      CGAL::Bbox_2 bb = c.bbox(); // assumes bounded curve
      int xmin = view->mapFromScene(bb.xmin(), bb.ymin()).x();
      int xmax = view->mapFromScene(bb.xmax(), bb.ymin()).x();
      n = xmax - xmin;
      if (n < 2) { n = 2; }
    }

    std::vector<std::pair<double, double>> app_pts;
    app_pts.reserve(n + 1);
    c.polyline_approximation(n, std::back_inserter(app_pts));
    auto end_pts = app_pts.end();
    auto p_curr = app_pts.begin();
    auto p_next = p_curr + 1;
    do
    {
      Point_2 p1(p_curr->first, p_curr->second);
      Point_2 p2(p_next->first, p_next->second);
      Segment_2 seg(p1, p2);

      FT dist = CGAL::squared_distance(p, seg);
      if (first || dist < min_dist)
      {
        first = false;
        min_dist = dist;
      }

      p_curr++;
      p_next++;
    } while (p_next != end_pts);

    return CGAL::to_double(min_dist);
  }
}

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
double Compute_squared_distance_2<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>>::
operator()(const Point_2& p, const X_monotone_curve_2& curve) const
{
  // TODO: this should probably be cached!
  CGAL::Qt::ArrangementPainterOstream<Traits> painterOstream{nullptr};
  painterOstream.setScene(this->getScene());

  std::pair<double, double> p_pair = {
    CGAL::to_double(p.x()), CGAL::to_double(p.y())};

  double minDist = (std::numeric_limits<double>::max)();
  auto points = painterOstream.getPoints(curve);
  for (auto& vit : points)
  {
    QPointF coord(vit.first, vit.second);
    float curDist = (vit.first - p_pair.first) * (vit.first - p_pair.first) +
                    (vit.second - p_pair.second) * (vit.second - p_pair.second);
    minDist = curDist < minDist ? curDist : minDist;
  }
  return minDist;
}

template <typename Coefficient_>
double
Compute_squared_distance_2<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
operator()(const Point_2& p, const X_monotone_curve_2& curve) const
{
  // TODO: this should probably be cached!
  CGAL::Qt::ArrangementPainterOstream<Traits> painterOstream{nullptr};
  painterOstream.setScene(this->getScene());

  auto view = this->getView();
  QTransform worldTransform;
  worldTransform.translate(
    -view->horizontalScrollBar()->value(), -view->verticalScrollBar()->value());
  worldTransform = view->transform() * worldTransform;

  QTransform facadeToViewport =
    painterOstream.getPointsListMapping(worldTransform);

  auto points = painterOstream.getPointsList(curve);

  QPoint p_viewport =
    view->mapFromScene(QPointF{p.x().doubleValue(), p.y().doubleValue()});

  double minDist = (std::numeric_limits<double>::max)();
  for (auto& vec : points)
  {
    for (auto vit = vec.begin(); vit != vec.end(); ++vit)
    {
      QPoint coord(vit->first, vit->second);
      float curDist = QLineF{facadeToViewport.map(coord), p_viewport}.length();
      minDist = curDist < minDist ? curDist : minDist;
    }
  }

  return minDist;
}

template <typename AlgebraicKernel_d_1>
double Compute_squared_distance_2<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>::
operator()(const Point_2& p, const X_monotone_curve_2& curve) const
{
  // TODO: this should probably be cached!
  CGAL::Qt::ArrangementPainterOstream<Traits> painterOstream{nullptr};
  painterOstream.setScene(this->getScene());

  std::pair<double, double> p_pair = {
    CGAL::to_double(p.x()), CGAL::to_double(p.y())};

  double minDist = (std::numeric_limits<double>::max)();
  auto points_list = painterOstream.getPointsList(curve);
  for (auto& points : points_list)
  {
    for (auto& vit : points)
    {
      QPointF coord(vit.first, vit.second);
      float curDist =
        (vit.first - p_pair.first) * (vit.first - p_pair.first) +
        (vit.second - p_pair.second) * (vit.second - p_pair.second);
      minDist = curDist < minDist ? curDist : minDist;
    }
  }
  return minDist;
}

template <typename ArrTraits>
auto Arr_construct_point_2<ArrTraits>::operator()(const Kernel_point_2& pt)
  -> Point_2
{
  return (*this)(FT{pt.x()}, FT{pt.y()}, traits);
}

template <typename ArrTraits>
auto Arr_construct_point_2<ArrTraits>::operator()(const FT& x, const FT& y)
  -> Point_2
{
  return (*this)(x, y, traits);
}

template <typename ArrTraits>
template <typename TTraits>
auto Arr_construct_point_2<ArrTraits>::operator()(
  const FT& x, const FT& y, const TTraits*) -> Point_2
{
  CoordinateType xx(x);
  CoordinateType yy(y);
  Point_2 res(xx, yy);
  return res;
}

template <typename ArrTraits>
template <typename AlgebraicKernel_d_1>
auto Arr_construct_point_2<ArrTraits>::operator()(
  const FT& x, const FT& y,
  const CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>*)
  -> Point_2
{
  using Rational = typename ArrTraits::Rational;
  using Rational_function = typename ArrTraits::Rational_function;
  using Polynomial_1 = typename ArrTraits::Polynomial_1;
  using RationalTraits = CGAL::Rational_traits<Rational>;

  static AlgebraicKernel_d_1 algebraic_kernel;

  Rational y_rat{y};
  RationalTraits ratTraits;

  CoordinateType xx(x);
  Polynomial_1 y_num{ratTraits.numerator(y_rat)};
  Polynomial_1 y_den{ratTraits.denominator(y_rat)};
  Rational_function yy{y_num, y_den, &algebraic_kernel};
  Point_2 res(yy, xx);
  return res;
}

template <typename Arr_>
auto Find_nearest_edge<Arr_>::operator()(const Point_2& queryPt)
  -> Halfedge_const_handle
{
  Face_const_handle face =
    PointLocationFunctions<Arrangement>{}.getFace(arr, queryPt);

  bool first = 1;
  X_monotone_curve_2 closestCurve;
  Halfedge_const_handle closestEdge;
  double minDist(0);

  if (!face->is_unbounded())
  { // it is an interior face so it has a ccb
    Ccb_halfedge_const_circulator cc = face->outer_ccb();
    do
    {
      X_monotone_curve_2 curve = cc->curve();
      double dist = this->pointCurveDistance(queryPt, curve);
      if (first || dist < minDist)
      {
        first = 0;
        minDist = dist;
        closestEdge = cc;
      }
    } while (++cc != face->outer_ccb());
  }
  else if (face->has_outer_ccb())
  {
    Ccb_halfedge_const_circulator cc = face->outer_ccb();
    do
    {
      if (cc->is_fictitious()) { continue; }

      X_monotone_curve_2 curve = cc->curve();
      double dist = this->pointCurveDistance(queryPt, curve);
      if (first || dist < minDist)
      {
        first = 0;
        minDist = dist;
        closestEdge = cc;
      }
    } while (++cc != face->outer_ccb());
  }
  Hole_const_iterator hit;
  Hole_const_iterator eit = face->holes_end();
  // int counter = 0;
  for (hit = face->holes_begin(); hit != eit; ++hit)
  { // check any holes inside this face
    Ccb_halfedge_const_circulator cc = *hit;
    do
    {
      X_monotone_curve_2 curve = cc->curve();
      double dist = this->pointCurveDistance(queryPt, curve);
      if (first || dist < minDist)
      {
        first = 0;
        minDist = dist;
        closestEdge = cc;
      }
      cc++;
    } while (cc != *hit);
  }

  return closestEdge;
}

template <typename Arr_>
auto Find_nearest_edge<Arr_>::getFace(const CGAL::Object& obj)
  -> Face_const_handle
{
  Face_const_handle f;
  if (CGAL::assign(f, obj)) return f;

  Halfedge_const_handle he;
  if (CGAL::assign(he, obj)) return (he->face());

  Vertex_const_handle v;
  CGAL_assertion(CGAL::assign(v, obj));
  CGAL::assign(v, obj);
  if (v->is_isolated()) return v->face();
  auto eit = v->incident_halfedges();
  return (eit->face());
}

template <typename ArrTraits>
Construct_x_monotone_subcurve_2<ArrTraits>::Construct_x_monotone_subcurve_2(
  const ArrTraits* traits_) :
    traits(traits_),
    split_2(this->traits->split_2_object()),
    compare_x_2(this->traits->compare_x_2_object()),
    compute_y_at_x(this->traits),
    construct_min_vertex_2(this->traits->construct_min_vertex_2_object()),
    construct_max_vertex_2(this->traits->construct_max_vertex_2_object()),
    parameter_space_in_x_2(this->traits)
{
}

template <typename ArrTraits>
auto Construct_x_monotone_subcurve_2<ArrTraits>::operator()(
  const X_monotone_curve_2& curve, const boost::optional<Point_2>& pLeft,
  const boost::optional<Point_2>& pRight) -> X_monotone_curve_2
{
  Point_2 pMin, pMax;
  bool unbounded_min = false;
  bool unbounded_max = false;

  if (this->parameter_space_in_x_2(curve, CGAL::ARR_MIN_END) != CGAL::INTERIOR)
    unbounded_min = true;
  else
    pMin = this->construct_min_vertex_2(curve);

  if (this->parameter_space_in_x_2(curve, CGAL::ARR_MAX_END) != CGAL::INTERIOR)
    unbounded_max = true;
  else
    pMax = this->construct_max_vertex_2(curve);

  X_monotone_curve_2 subcurve;
  X_monotone_curve_2 unusedTrimmings;
  X_monotone_curve_2 finalSubcurve;
  if (
    pLeft && (unbounded_min || this->compare_x_2(*pLeft, pMin) == CGAL::LARGER))
  {
    auto y1 = this->compute_y_at_x(curve, pLeft->x());

    Point_2 splitPoint = {pLeft->x(), y1};
    this->split_2(curve, splitPoint, unusedTrimmings, subcurve);
  }
  else
  {
    subcurve = curve;
  }

  if (
    pRight &&
    (unbounded_max || this->compare_x_2(*pRight, pMax) == CGAL::SMALLER))
  {
    auto y2 = this->compute_y_at_x(subcurve, pRight->x());
    Point_2 splitPoint = {pRight->x(), y2};
    this->split_2(subcurve, splitPoint, finalSubcurve, unusedTrimmings);
  }
  else
  {
    finalSubcurve = subcurve;
  }

  return finalSubcurve;
}

template <typename RatKernel, typename AlgKernel, typename NtTraits>
auto Construct_x_monotone_subcurve_2<
  CGAL::Arr_conic_traits_2<RatKernel, AlgKernel, NtTraits>>::
operator()(
  const X_monotone_curve_2& curve, const boost::optional<Point_2>& pLeft,
  const boost::optional<Point_2>& pRight) -> X_monotone_curve_2
{
  // TODO: handle when pLeft or pRight is null

  // find the points on the curve
  Point_2 left = curve.point_at_x(*pLeft);
  Point_2 right = curve.point_at_x(*pRight);

  // make sure the points are oriented in the direction that the curve is
  // going
  AlgKernel ker;
  if (!(((curve.is_directed_right()) &&
         ker.compare_xy_2_object()(left, right) == CGAL::SMALLER) ||
        ((!curve.is_directed_right()) &&
         ker.compare_xy_2_object()(left, right) == CGAL::LARGER)))
  { std::swap(left, right); }

  X_monotone_curve_2 res = curve.trim(left, right);
  return res;
}

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
Construct_x_monotone_subcurve_2<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits,
  BoundingTraits>>::Construct_x_monotone_subcurve_2(const ArrTraits* traits_) :
    traits(traits_),
    split_2(this->traits->split_2_object()),
    compare_x_2(this->traits->compare_x_2_object()),
    compute_y_at_x(this->traits),
    construct_min_vertex_2(this->traits->construct_min_vertex_2_object()),
    construct_max_vertex_2(this->traits->construct_max_vertex_2_object())
{
}

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
auto Construct_x_monotone_subcurve_2<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>>::
operator()(
  const X_monotone_curve_2& curve, const boost::optional<Point_2>& pLeft,
  const boost::optional<Point_2>& pRight) -> X_monotone_curve_2
{
  auto pMin = this->construct_min_vertex_2(curve);
  auto pMax = this->construct_max_vertex_2(curve);

  X_monotone_curve_2 subcurve;
  X_monotone_curve_2 unusedTrimmings;
  X_monotone_curve_2 finalSubcurve;

  // FIXME (maybe no way to fix it?): This is not optimal
  // the envelope package returns x-monotone curves, and range [xbegin, xend]
  // the range is algebraic, but we can't find "t" of the bezier curve
  // at algebraic "x", so we need to make it rational
  auto local_get_t = [&](auto&& point) {
    if (point.is_rational())
      return this->compute_y_at_x.get_t(
        curve, ((typename Point_2::Rat_point_2)point).x());
    else
      return this->compute_y_at_x.get_t(curve, point.approximate().first);
  };

  if (pLeft && this->compare_x_2(*pLeft, pMin) == CGAL::LARGER)
  {
    auto t = local_get_t(*pLeft);
    Point_2 splitPoint(curve.supporting_curve(), t);
    this->split_2(curve, splitPoint, unusedTrimmings, subcurve);
  }
  else
  {
    subcurve = curve;
  }
  if (pRight && this->compare_x_2(*pRight, pMax) == CGAL::SMALLER)
  {
    auto t = local_get_t(*pRight);
    Point_2 splitPoint(curve.supporting_curve(), t);
    this->split_2(subcurve, splitPoint, finalSubcurve, unusedTrimmings);
  }
  else
  {
    finalSubcurve = subcurve;
  }

  return finalSubcurve;
}

template <typename AlgebraicKernel_d_1>
Construct_x_monotone_subcurve_2<CGAL::Arr_rational_function_traits_2<
  AlgebraicKernel_d_1>>::Construct_x_monotone_subcurve_2(const Traits*
                                                           traits_) :
    traits(traits_),
    split_2(this->traits->split_2_object()),
    compare_x_2(this->traits->compare_x_2_object()),
    compute_y_at_x(this->traits),
    construct_min_vertex_2(this->traits->construct_min_vertex_2_object()),
    construct_max_vertex_2(this->traits->construct_max_vertex_2_object())
{
}

template <typename AlgebraicKernel_d_1>
auto Construct_x_monotone_subcurve_2<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>::
operator()(
  const X_monotone_curve_2& curve, const boost::optional<Point_2>& pLeft,
  const boost::optional<Point_2>& pRight) -> X_monotone_curve_2
{
  Point_2 pMin, pMax;
  bool unbounded_min = false;
  bool unbounded_max = false;
  try { pMin = this->construct_min_vertex_2(curve); }
  catch (...) { unbounded_min = true; }
  try { pMax = this->construct_max_vertex_2(curve); }
  catch (...) { unbounded_max = true; }

  X_monotone_curve_2 subcurve;
  X_monotone_curve_2 unusedTrimmings;
  X_monotone_curve_2 finalSubcurve;
  if (
    pLeft && (unbounded_min || this->compare_x_2(*pLeft, pMin) == CGAL::LARGER))
  {
    Point_2 splitPoint{curve._f, pLeft->x()};
    this->split_2(curve, splitPoint, unusedTrimmings, subcurve);
  }
  else
  {
    subcurve = curve;
  }

  if (
    pRight &&
    (unbounded_max || this->compare_x_2(*pRight, pMax) == CGAL::SMALLER))
  {
    auto y2 = this->compute_y_at_x(subcurve, pRight->x());
    Point_2 splitPoint{curve._f, pRight->x()};
    this->split_2(subcurve, splitPoint, finalSubcurve, unusedTrimmings);
  }
  else
  {
    finalSubcurve = subcurve;
  }

  return finalSubcurve;
}

template <typename ArrTraits>
Arr_compute_y_at_x_2<ArrTraits>::Arr_compute_y_at_x_2(const Traits* traits_) :
    traits(traits_),
    intersectCurves(this->traits->intersect_2_object())
{
}

template <typename ArrTraits>
auto Arr_compute_y_at_x_2<ArrTraits>::operator()(
  const X_monotone_curve_2& curve, const CoordinateType& x) -> CoordinateType
{
  typename Traits::Left_side_category category;
  return this->operator()(curve, x, this->traits, category);
}

template <typename ArrTraits>
double Arr_compute_y_at_x_2<ArrTraits>::approx(
  const X_monotone_curve_2& curve, const CoordinateType& x)
{
  return CGAL::to_double((*this)(curve, x));
}

template <typename ArrTraits>
template <typename TTraits>
auto Arr_compute_y_at_x_2<ArrTraits>::operator()(
  const X_monotone_curve_2& curve, const CoordinateType& x,
  const TTraits* traits_, CGAL::Arr_oblivious_side_tag) -> CoordinateType
{
  typedef
    typename TTraits::Construct_x_monotone_curve_2 Construct_x_monotone_curve_2;
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2 =
    traits_->construct_x_monotone_curve_2_object();
  CoordinateType res(0);
  CGAL::Bbox_2 clipRect = curve.bbox();
  Point_2 p1c1(x, CoordinateType(clipRect.ymin() - 1)); // clicked point
  // upper bounding box
  Point_2 p2c1(x, CoordinateType(clipRect.ymax() + 1));

  const X_monotone_curve_2 verticalLine =
    construct_x_monotone_curve_2(p1c1, p2c1);
  CGAL::Object o;
  CGAL::Oneset_iterator<CGAL::Object> oi(o);

  this->intersectCurves(curve, verticalLine, oi);

  IntersectionResult pair;
  if (CGAL::assign(pair, o))
  {
    Point_2 pt = pair.first;
    res = pt.y();
  }
  return res;
}

template <typename ArrTraits>
template <typename TTraits>
auto Arr_compute_y_at_x_2<ArrTraits>::operator()(
  const X_monotone_curve_2& curve, const CoordinateType& x,
  const TTraits* traits_, CGAL::Arr_open_side_tag) -> CoordinateType
{
  typename TTraits::Construct_x_monotone_curve_2 construct_x_monotone_curve_2 =
    traits_->construct_x_monotone_curve_2_object();
  CoordinateType res(0);
  // QRectF clipRect = this->viewportRect( );
  Line_2 line = curve.supporting_line();
  // FIXME: get a better bounding box for an unbounded segment
  Point_2 p1c1(x, CoordinateType(-10000000)); // clicked point
  Point_2 p2c1(x, CoordinateType(10000000));  // upper bounding box

  const X_monotone_curve_2 verticalLine =
    construct_x_monotone_curve_2(p1c1, p2c1);
  CGAL::Object o;
  CGAL::Oneset_iterator<CGAL::Object> oi(o);

  this->intersectCurves(curve, verticalLine, oi);

  IntersectionResult pair;
  if (CGAL::assign(pair, o))
  {
    Point_2 pt = pair.first;
    res = pt.y();
  }
  return res;
}

template <typename Coefficient_>
auto Arr_compute_y_at_x_2<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
operator()(const X_monotone_curve_2& curve, const CoordinateType& x)
  -> CoordinateType
{
  CGAL::Object o;
  CGAL::Oneset_iterator<CGAL::Object> oi(o);
  Intersect_2 intersect = traits->intersect_2_object();
  X_monotone_curve_2 c2 = this->makeVerticalLine(x);
  intersect(curve, c2, oi);
  std::pair<Point_2, Multiplicity> res;
  if (CGAL::assign(res, o)) // TODO: handle failure case
  {
    const Point_2& p = res.first;
    CoordinateType coord = p.y();
    return coord;
  }
  else
  {
    std::cout << "Warning: vertical projection failed" << std::endl;
    return CoordinateType(0);
  }
}

template <typename Coefficient_>
double
Arr_compute_y_at_x_2<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  approx(const X_monotone_curve_2& curve, const CoordinateType& x)
{
  return CGAL::to_double(this->operator()(curve, x));
}

template <typename Coefficient_>
auto Arr_compute_y_at_x_2<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
  makeVerticalLine(const CoordinateType& x) -> X_monotone_curve_2
{
  typename Traits::Construct_point_2 constructPoint =
    traits->construct_point_2_object();
  typename Traits::Construct_x_monotone_segment_2 constructSegment =
    traits->construct_x_monotone_segment_2_object();

  std::vector<X_monotone_curve_2> curves;
  Point_2 p1 =
    constructPoint(x, CoordinateType(-(std::numeric_limits<int>::max)()));
  Point_2 p2 =
    constructPoint(x, CoordinateType((std::numeric_limits<int>::max)()));
  constructSegment(p1, p2, std::back_inserter(curves));
  return curves[0]; // by construction, there is one curve in curves
}

template <typename Bezier_x_monotone_2>
static inline auto get_t_range(const Bezier_x_monotone_2& curve)
{
  auto&& supp_curve = curve.supporting_curve();
  auto ps_org = curve.source().get_originator(supp_curve, curve.xid());
  CGAL_assertion(ps_org != curve.source().originators_end());

  auto pt_org = curve.target().get_originator(supp_curve, curve.xid());
  CGAL_assertion(pt_org != curve.target().originators_end());

  return std::make_pair(
    (ps_org->point_bound().t_min + ps_org->point_bound().t_max) / 2,
    (pt_org->point_bound().t_min + pt_org->point_bound().t_max) / 2);
}

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
auto Arr_compute_y_at_x_2<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>>::
operator()(const X_monotone_curve_2& curve, const Rational& x) -> Algebraic
{
  auto&& supp_curve = curve.supporting_curve();
  auto t = this->get_t(curve, x);
  NtTraits nt_traits;
  return nt_traits.evaluate_at(supp_curve.y_polynomial(), t) /
         nt_traits.convert(supp_curve.y_norm());
}

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
auto Arr_compute_y_at_x_2<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits, BoundingTraits>>::
  get_t(const X_monotone_curve_2& curve, const Rational& x) -> Algebraic
{
  std::vector<Algebraic> t_vals;
  curve.supporting_curve().get_t_at_x(x, std::back_inserter(t_vals));
  auto t_range = get_t_range(curve);

  const Algebraic& t_src{t_range.first};
  const Algebraic& t_trg{t_range.second};

  for (auto t_iter = t_vals.begin(); t_iter != t_vals.end(); ++t_iter)
  {
    auto res1 = CGAL::compare(*t_iter, t_src);
    if (res1 == CGAL::EQUAL) return (curve.source().y());
    auto res2 = CGAL::compare(*t_iter, t_trg);
    if (res2 == CGAL::EQUAL) return (curve.target().y());
    if (res1 != res2) return *t_iter;
  }

  CGAL_error();
  return 0;
}

template <
  typename RatKernel, typename AlgKernel, typename NtTraits,
  typename BoundingTraits>
double Arr_compute_y_at_x_2<CGAL::Arr_Bezier_curve_traits_2<
  RatKernel, AlgKernel, NtTraits,
  BoundingTraits>>::approx(const X_monotone_curve_2& curve, const Rational& x)
{
  return CGAL::to_double((*this)(curve, x));
}

template <typename AlgebraicKernel_d_1>
auto Arr_compute_y_at_x_2<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>::
operator()(const X_monotone_curve_2& curve, const Algebraic_real_1& x)
  -> Algebraic_real_1
{
  return Point_2{curve._f, x}.y();
}

template <typename AlgebraicKernel_d_1>
auto Arr_compute_y_at_x_2<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>::
operator()(const X_monotone_curve_2& curve, const Rational& x) -> Rational
{
  return curve._f.numer().evaluate(x) / curve._f.denom().evaluate(x);
}

template <typename AlgebraicKernel_d_1>
auto Arr_compute_y_at_x_2<
  CGAL::Arr_rational_function_traits_2<AlgebraicKernel_d_1>>::
approx(const X_monotone_curve_2& curve, const Rational& x) -> double
{
  return CGAL::to_double((*this)(curve, x));
}

CGAL::Object createArrangement(demo_types::TraitsType tt)
{
  CGAL::Object res;
  demo_types::visitArrangementType(tt, [&](auto type_holder) {
    using Arrangement = typename decltype(type_holder)::type;
    res = CGAL::make_object(new Arrangement());
  });
  return res;
}

// Insert_curve
template <typename Arr_>
void Insert_curve<Arr_>::operator()(Arrangement* arr, const Curve_2& curve)
{
    CGAL::insert(*arr, curve);
}

// free functions
void deleteArrangement(demo_types::TraitsType tt, const CGAL::Object& arr_obj)
{
  demo_types::visitArrangementType(tt, [&](auto type_holder) {
    using Arrangement = typename decltype(type_holder)::type;
    Arrangement* arr = nullptr;
    CGAL::assign(arr, arr_obj);
    delete arr;
  });
}

CGAL::Object makeOverlayArrangement(const std::vector<CGAL::Object>& arrs)
{
  CGAL::Object arr_obj;
  if (arrs.size() == 2)
  {
    demo_types::forEachArrangementType([&](auto type_holder) {
      using Arrangement = typename decltype(type_holder)::type;

      Arrangement* arr1;
      Arrangement* arr2;
      if (CGAL::assign(arr1, arrs[0]) && CGAL::assign(arr2, arrs[1]))
      {
        auto overlay_arr = new Arrangement();
        CGAL::Arr_default_overlay_traits<Arrangement> overlay_traits;

        CGAL::overlay(*arr1, *arr2, *overlay_arr, overlay_traits);
        arr_obj = CGAL::make_object(overlay_arr);
      }
    });
  }
  return arr_obj;
}

void insertCurve(
  demo_types::TraitsType tt, const CGAL::Object& arr_obj,
  const CGAL::Object& curve_obj)
{
  demo_types::visitArrangementType(tt, [&](auto type_holder) {
    using Arrangement = typename decltype(type_holder)::type;
    using Curve_2 = typename Arrangement::Curve_2;

    Curve_2 curve;
    Arrangement* arr;
    if (!CGAL::assign(arr, arr_obj)) CGAL_error();
    if (!CGAL::assign(curve, curve_obj)) CGAL_error();

    Insert_curve<Arrangement>{}(arr, curve);
  });
}

ARRANGEMENT_DEMO_SPECIALIZE_TRAITS(Compute_squared_distance_2)
ARRANGEMENT_DEMO_SPECIALIZE_TRAITS(Arr_construct_point_2)
ARRANGEMENT_DEMO_SPECIALIZE_TRAITS(Construct_x_monotone_subcurve_2)
ARRANGEMENT_DEMO_SPECIALIZE_TRAITS(Arr_compute_y_at_x_2)
ARRANGEMENT_DEMO_SPECIALIZE_ARR(Find_nearest_edge)
ARRANGEMENT_DEMO_SPECIALIZE_ARR(Insert_curve)
