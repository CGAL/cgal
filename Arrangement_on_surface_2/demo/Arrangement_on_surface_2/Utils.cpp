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

#include "Utils.h"
#include "ArrangementTypes.h"
#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>

#include <QScrollBar>

template <typename Kernel_>
double
Compute_squared_distance_2<CGAL::Arr_segment_traits_2<Kernel_>>::operator()(
  const Point_2& p, const X_monotone_curve_2& c) const
{
  Point_2 p1 = c.source();
  Point_2 p2 = c.target();
  Segment_2 seg(p1, p2);

  return CGAL::to_double(this->squared_distance(p, seg));
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
    res = this->squared_distance(p, seg);
  }
  else if (c.is_ray())
  {
    ray = c.ray();
    res = this->squared_distance(p, ray);
  }
  else // ( c.is_line( ) )
  {
    line = c.line();
    res = this->squared_distance(p, line);
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
    FT dist = this->squared_distance(p, seg);

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

    std::pair<double, double>* app_pts = new std::pair<double, double>[n + 1];
    std::pair<double, double>* end_pts = c.polyline_approximation(n, app_pts);
    std::pair<double, double>* p_curr = app_pts;
    std::pair<double, double>* p_next = p_curr + 1;
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

template <typename RatKernel, typename AlgKernel, typename NtTraits>
double Compute_squared_distance_2<
  CGAL::Arr_Bezier_curve_traits_2<RatKernel, AlgKernel, NtTraits>>::
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

    std::pair<double, double>* app_pts = new std::pair<double, double>[n + 1];
    std::pair<double, double>* end_pts = c.polyline_approximation(n, app_pts);
    std::pair<double, double>* p_curr = app_pts;
    std::pair<double, double>* p_next = p_curr + 1;
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

template <typename Coefficient_>
double
Compute_squared_distance_2<CGAL::Arr_algebraic_segment_traits_2<Coefficient_>>::
operator()(const Point_2& p, const X_monotone_curve_2& c) const
{
  // TODO: this should probably be cached!
  // this is the same as ArrangementPainterOstream::setupFacade
  // TODO: refactor it
  typedef CGAL::Curve_renderer_facade<CKvA_2> Facade;
  QGraphicsView* view = this->getView();
  QRectF viewport = this->viewportRect();
  CGAL::Bbox_2 bbox = CGAL::Qt::Converter<Kernel>{}(viewport).bbox();
  Facade::setup(bbox, view->width(), view->height());

  // this is the same as ArrangementPainterOstream::remapFacadePainter
  // TODO: refactor it

  // this is equivalent to QPainter::worldTransform
  QTransform worldTransform;
  worldTransform.translate(
    -view->horizontalScrollBar()->value(), -view->verticalScrollBar()->value());
  worldTransform = view->transform() * worldTransform;

  QPointF dxdy = worldTransform.map(viewport.topLeft());
  QPointF p1 = worldTransform.map(viewport.topRight());
  QPointF p2 = worldTransform.map(viewport.bottomLeft());
  float dx = dxdy.x();
  float dy = dxdy.y();
  float m11 = (p1.x() - dx) / view->width();
  float m21 = (p2.x() - dx) / view->height();
  float m22 = (p2.y() - dy) / view->height();
  float m12 = (p1.y() - dy) / view->width();
  auto facadeToViewport = QTransform{m11, m12, m21, m22, dx, dy};

  std::list<Coord_vec_2> points;
  Facade::instance().draw(c, points);

  QPoint p_viewport =
    view->mapFromScene(QPointF{p.x().doubleValue(), p.y().doubleValue()});

  double minDist = std::numeric_limits<double>::max();
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

template <typename ArrTraits>
auto Arr_construct_point_2<ArrTraits>::operator()(const Kernel_point_2& pt)
  -> Point_2
{
  return (*this)(pt.x(), pt.y());
}

template <typename ArrTraits>
template <typename T>
auto Arr_construct_point_2<ArrTraits>::operator()(const T& x, const T& y)
  -> Point_2
{
  return (*this)(x, y, ArrTraits());
}

template <typename ArrTraits>
template <typename T, typename TTraits>
auto Arr_construct_point_2<ArrTraits>::operator()(
  const T& x, const T& y, TTraits /* traits */) -> Point_2
{
  CoordinateType xx(x);
  CoordinateType yy(y);
  Point_2 res(xx, yy);
  return res;
}

template <typename Arr_, typename ArrTraits>
auto Find_nearest_edge<Arr_, ArrTraits>::operator()(const Point_2& queryPt)
  -> Halfedge_const_handle
{
  typename ArrTraits::Point_2 pt = this->toArrPoint(queryPt);
  CGAL::Object pointLocationResult = this->pointLocationStrategy.locate(pt);
  Face_const_handle face = this->getFace(pointLocationResult);
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

template <typename Arr_, typename ArrTraits>
auto Find_nearest_edge<Arr_, ArrTraits>::getFace(const CGAL::Object& obj)
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
  Halfedge_around_vertex_const_circulator eit = v->incident_halfedges();
  return (eit->face());
}

template <typename ArrTraits>
Construct_x_monotone_subcurve_2<ArrTraits>::Construct_x_monotone_subcurve_2() :
    intersect_2(this->traits.intersect_2_object()),
    split_2(this->traits.split_2_object()),
    compare_x_2(this->traits.compare_x_2_object()),
    construct_min_vertex_2(this->traits.construct_min_vertex_2_object()),
    construct_max_vertex_2(this->traits.construct_max_vertex_2_object())
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
  try
  {
    pMin = this->construct_min_vertex_2(curve);
  }
  catch (...)
  {
    unbounded_min = true;
  }
  try
  {
    pMax = this->construct_max_vertex_2(curve);
  }
  catch (...)
  {
    unbounded_max = true;
  }

  X_monotone_curve_2 subcurve;
  X_monotone_curve_2 unusedTrimmings;
  X_monotone_curve_2 finalSubcurve;
  if (
    pLeft && (unbounded_min || this->compare_x_2(*pLeft, pMin) == CGAL::LARGER))
  {
    CoordinateType y1 = this->compute_y_at_x(curve, pLeft->x());
    Point_2 splitPoint(pLeft->x(), y1);
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
    CoordinateType y2 = this->compute_y_at_x(subcurve, pRight->x());
    Point_2 splitPoint(pRight->x(), y2);
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

template <typename ArrTraits>
Arr_compute_y_at_x_2<ArrTraits>::Arr_compute_y_at_x_2() :
    intersectCurves(this->traits.intersect_2_object())
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
  const X_monotone_curve_2& curve, const CoordinateType& x, TTraits traits_,
  CGAL::Arr_oblivious_side_tag) -> CoordinateType
{
  typedef
    typename TTraits::Construct_x_monotone_curve_2 Construct_x_monotone_curve_2;
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2 =
    traits_.construct_x_monotone_curve_2_object();
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
  const X_monotone_curve_2& curve, const CoordinateType& x, TTraits traits_,
  CGAL::Arr_open_side_tag) -> CoordinateType
{
  typename TTraits::Construct_x_monotone_curve_2 construct_x_monotone_curve_2 =
    traits_.construct_x_monotone_curve_2_object();
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
operator()(
  const X_monotone_curve_2& curve, const CoordinateType& x, Point_2* out)
  -> CoordinateType
{
  CGAL::Object o;
  CGAL::Oneset_iterator<CGAL::Object> oi(o);
  Intersect_2 intersect = traits.intersect_2_object();
  X_monotone_curve_2 c2 = this->makeVerticalLine(x);
  intersect(curve, c2, oi);
  std::pair<Point_2, Multiplicity> res;
  if (CGAL::assign(res, o)) // TODO: handle failure case
  {
    Point_2 p = res.first;
    CoordinateType coord = p.y();
    if (out) *out = p;
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
    traits.construct_point_2_object();
  typename Traits::Construct_x_monotone_segment_2 constructSegment =
    traits.construct_x_monotone_segment_2_object();

  std::vector<X_monotone_curve_2> curves;
  Point_2 p1 = constructPoint(x, CoordinateType(-1000000));
  Point_2 p2 = constructPoint(x, CoordinateType(+1000000));
  constructSegment(p1, p2, std::back_inserter(curves));
  return curves[0]; // by construction, there is one curve in curves
}

// TODO(Ahmed Essam): create a macro to specialize structs
template class Compute_squared_distance_2<Seg_traits>;
template class Compute_squared_distance_2<Pol_traits>;
template class Compute_squared_distance_2<Conic_traits>;
template class Compute_squared_distance_2<Lin_traits>;
template class Compute_squared_distance_2<Alg_seg_traits>;
// template class Compute_squared_distance_2<Bezier_traits>;

template class Arr_construct_point_2<Seg_traits>;
template class Arr_construct_point_2<Pol_traits>;
template class Arr_construct_point_2<Conic_traits>;
template class Arr_construct_point_2<Lin_traits>;
template class Arr_construct_point_2<Alg_seg_traits>;
// template class Arr_construct_point_2<Bezier_traits>;

template class Find_nearest_edge<Seg_arr, Seg_traits>;
template class Find_nearest_edge<Pol_arr, Pol_traits>;
template class Find_nearest_edge<Conic_arr, Conic_traits>;
template class Find_nearest_edge<Lin_arr, Lin_traits>;
template class Find_nearest_edge<Alg_seg_arr, Alg_seg_traits>;
// template class Find_nearest_edge<Bezier_arr, Bezier_traits>;

template class Construct_x_monotone_subcurve_2<Seg_traits>;
template class Construct_x_monotone_subcurve_2<Pol_traits>;
template class Construct_x_monotone_subcurve_2<Conic_traits>;
template class Construct_x_monotone_subcurve_2<Lin_traits>;
template class Construct_x_monotone_subcurve_2<Alg_seg_traits>;
// template class Construct_x_monotone_subcurve_2<Bezier_traits>;

template class Arr_compute_y_at_x_2<Seg_traits>;
template class Arr_compute_y_at_x_2<Pol_traits>;
template class Arr_compute_y_at_x_2<Conic_traits>;
template class Arr_compute_y_at_x_2<Lin_traits>;
template class Arr_compute_y_at_x_2<Alg_seg_traits>;
// template class Arr_compute_y_at_x_2<Bezier_traits>;
