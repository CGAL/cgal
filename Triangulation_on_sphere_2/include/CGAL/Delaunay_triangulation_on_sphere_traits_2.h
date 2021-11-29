// Copyright (c) 1997, 2012, 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mariette Yvinec,
//                 Claudia Werner,
//                 Mael Rouxel-Labbé

#ifndef CGAL_DELAUNAY_TRIANGULATION_SPHERE_TRAITS_2_H
#define CGAL_DELAUNAY_TRIANGULATION_SPHERE_TRAITS_2_H

#include <CGAL/license/Triangulation_on_sphere_2.h>

#include <CGAL/Triangulation_on_sphere_2/internal/get_precision_bounds.h>

#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Has_conversion.h>
#include <CGAL/Spherical_kernel_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/triangulation_assertions.h>

#include <CGAL/enum.h>

namespace CGAL {
namespace internal {

template <typename LK>
class Orientation_on_sphere_2
{
public:
  typedef typename LK::Point_3                      Point_3;
  typedef typename LK::Point_3                      Point_on_sphere_2;
  typedef Comparison_result                         result_type;

  Orientation_on_sphere_2(const Point_3& center, const LK& lk)
    : _center(center), _lk(lk)
  { }

  Comparison_result operator()(const Point_on_sphere_2& p, const Point_on_sphere_2& q, const Point_on_sphere_2& r) const
  { return _lk.orientation_3_object()(_center, p, q, r); }

protected:
  const Point_3& _center;
  const LK& _lk;
};

template <typename LK>
class Equal_on_sphere_2
{
public:
  typedef typename LK::Point_3                      Point_3;
  typedef typename LK::Point_3                      Point_on_sphere_2;
  typedef bool                                      result_type;

  Equal_on_sphere_2(const Point_3& center, const LK& lk)
    : _center(center), _lk(lk)
  { }

  bool operator()(const Point_on_sphere_2& p, const Point_on_sphere_2& q) const
  {
    return _lk.collinear_3_object()(_center, p, q) &&
            !_lk.collinear_are_strictly_ordered_along_line_3_object()(p, _center, q);
  }

protected:
  const Point_3& _center;
  const LK& _lk;
};

template <typename LK>
class Inside_cone_2
{
public:
  typedef typename LK::Point_3                      Point_3;
  typedef typename LK::Point_3                      Point_on_sphere_2;
  typedef bool                                      result_type;

  Inside_cone_2(const Point_3& center, const LK& lk)
    : _center(center), _lk(lk)
  { }

  bool operator()(const Point_on_sphere_2& p, const Point_on_sphere_2& q, const Point_on_sphere_2& r) const
  {
    // @todo first two checks might be unnecessary because we should always have r != p|q (using equal_on_sphere)
    if(_lk.collinear_3_object()(_center, p, r) ||
       _lk.collinear_3_object()(_center, q, r) ||
       _lk.orientation_3_object()(_center, p, q, r) != COLLINEAR)
    {
      return false;
    }

    // @todo also useless?
    if(collinear(_center, p, q))
      return true;

    // Check if r is on the same side of the plane orthogonal to the plane (c,p,q)
    // and going through c and p as q
    const Orientation op = _lk.coplanar_orientation_3_object()(_center, p, q, r);
    CGAL_triangulation_assertion(op != COLLINEAR);
    if(op == NEGATIVE)
      return false;

    // Check if r is on the same side of the plane orthogonal to the plane (c,p,q)
    // and going through c and q as p
    const Orientation oq = _lk.coplanar_orientation_3_object()(_center, q, p, r);
    CGAL_triangulation_assertion(oq != COLLINEAR);
    if(oq == NEGATIVE)
      return false;

    return true; // both true
  }

protected:
  const Point_3& _center;
  const LK& _lk;
};

template <typename LK>
typename LK::Point_3
ray_sphere_intersection(const typename LK::Point_3 o, // origin of the ray
                        const typename LK::Vector_3 dir, // direction of the ray
                        const typename LK::Point_3 sc, // center of the sphere
                        const typename LK::FT r, // radius of the sphere
                        const LK& lk)
{
  typedef typename LK::FT                                            FT;
  typedef typename LK::Vector_3                                      Vector_3;

  // Line-Sphere intersection is:
  //
  //   || x - c ||^2 = r^2    // sphere
  //   x = o + t * d          // ray
  //
  // thus || o + t*d - c ||^2 = r^2
  //
  // so a*t^2 + b*t + c = 0 with
  //
  // a := || d ||^2
  // b := 2 * d * (o - c)
  // c := || o - c||^2 - r^2

  const Vector_3 sco = lk.construct_vector_3_object()(sc, o);
  const FT a = lk.compute_squared_length_3_object()(dir);
  const FT b = 2 * lk.compute_scalar_product_3_object()(dir, sco);
  const FT c = lk.compute_squared_length_3_object()(sco) - r*r;

  FT d = b*b - 4*a*c;
  // by construction d > 0 is expected, but just to be safe
  if(d < 0)
    d = 0;

  const FT dr = CGAL::sqrt(d);
  const FT det = FT(0.5) / a;

  // care only about the positive value, and we know det > 0 && dr > 0
  const FT t0 = det * (-b + dr);

  return lk.construct_translated_point_3_object()(o, lk.construct_scaled_vector_3_object()(dir, t0));
}

template <typename LK>
class Construct_circumcenter_on_sphere_2
{
public:
  typedef typename LK::FT                       FT;
  typedef typename LK::Point_3                  Point_3;
  typedef typename LK::Point_3                  Point_on_sphere_2;
  typedef typename LK::Vector_3                 Vector_3;
  typedef Point_on_sphere_2                     result_type;

  Construct_circumcenter_on_sphere_2(const Point_3& center, const FT radius, const LK& lk)
    : _center(center), _radius(radius), _lk(lk)
  { }

  result_type operator()(const Point_on_sphere_2& p, const Point_on_sphere_2& q, const Point_on_sphere_2& r) const
  {
    const Point_3 c = _lk.construct_circumcenter_3_object()(p, q, r);
    const Vector_3 n = _lk.construct_normal_3_object()(p, q, r);
    return ray_sphere_intersection(c, n, _center, _radius, _lk);
  }

protected:
  const Point_3& _center;
  const FT _radius;
  const LK& _lk;
};

template <typename LK, typename SK>
class Construct_arc_on_sphere_2
{
public:
  typedef typename LK::FT                       FT;
  typedef typename LK::Point_3                  Point_3;
  typedef typename LK::Point_3                  Point_on_sphere_2;
  typedef typename SK::Circular_arc_3           result_type;

  Construct_arc_on_sphere_2(const Point_3& center, const FT radius, const LK& lk, const SK& sk)
    : _center(center), _radius(radius), _lk(lk), _sk(sk)
  { }

  result_type operator()(const Point_on_sphere_2& p, const Point_on_sphere_2& q) const
  {
    typedef typename CGAL::internal::Converter_selector<LK, SK>::type Converter;
    Converter cv;

    typename SK::Point_3 sc = cv(_center); // @tmp cache that
    typename SK::Point_3 sp = cv(p);
    typename SK::Point_3 sq = cv(q);

    CGAL_precondition(!_sk.collinear_3_object()(sc, sp, sq));

    typename SK::Plane_3 pl = _sk.construct_plane_3_object()(sc, sp, sq);
    CGAL_assertion(is_valid(pl));

    // CircleC3 has an annoying assertion to check pl.has_on(sc) that has no chance to
    // pass with a non-exact kernel, so dodging it by passing the plane's normal,
    //despite this creating internally another plane...
    typename SK::Circle_3 c = _sk.construct_circle_3_object()(
                                sc, square(_radius), _sk.construct_orthogonal_vector_3_object()(pl));

    typename SK::Circular_arc_point_3 cp = _sk.construct_circular_arc_point_3_object()(sp);
    typename SK::Circular_arc_point_3 cq = _sk.construct_circular_arc_point_3_object()(sq);

    return _sk.construct_circular_arc_3_object()(c, cp, cq);
  }

protected:
  const Point_3& _center;
  const FT _radius;
  const LK& _lk;
  const SK& _sk;
};

} // namespace internal

template <typename LK_,
          typename SK_ = CGAL::Spherical_kernel_3<
                          LK_, CGAL::Algebraic_kernel_for_spheres_2_3<typename LK_::FT> > >
class Delaunay_triangulation_on_sphere_traits_2
  : public LK_, SK_
{
public:
  typedef LK_                                                        LK;
  typedef SK_                                                        SK;

private:
  typedef Delaunay_triangulation_on_sphere_traits_2<LK>              Self;

public:
  typedef typename LK::FT                                            FT;
  typedef typename LK::Point_3                                       Point_on_sphere_2;
  typedef typename LK::Point_3                                       Point_3;
  typedef typename LK::Segment_3                                     Segment_3;
  typedef typename LK::Triangle_3                                    Triangle_3;

  typedef typename SK::Circular_arc_3                                Arc_on_sphere_2;

  // predicates
  typedef internal::Inside_cone_2<LK>                                Collinear_are_strictly_ordered_on_great_circle_2;
  typedef typename LK::Compare_xyz_3                                 Compare_on_sphere_2;
  typedef internal::Equal_on_sphere_2<LK>                            Equal_on_sphere_2;
  typedef internal::Orientation_on_sphere_2<LK>                      Orientation_on_sphere_2;
  typedef typename LK::Orientation_3                                 Side_of_oriented_circle_on_sphere_2;

  // constructions
  typedef typename LK::Construct_point_3                             Construct_point_on_sphere_2;
  typedef typename LK::Construct_point_3                             Construct_point_3;
  typedef typename LK::Construct_segment_3                           Construct_segment_3;
  typedef typename LK::Construct_triangle_3                          Construct_triangle_3;

  // For the hilbert sort
  typedef typename LK::Compute_x_3                                   Compute_x_3;
  typedef typename LK::Compute_y_3                                   Compute_y_3;
  typedef typename LK::Compute_z_3                                   Compute_z_3;

  typedef internal::Construct_circumcenter_on_sphere_2<LK>           Construct_circumcenter_on_sphere_2;
  typedef internal::Construct_arc_on_sphere_2<LK, SK>                Construct_arc_on_sphere_2;
  typedef typename LK::Construct_circumcenter_3                      Construct_circumcenter_3;

protected:
  Point_3 _center;
  FT _radius;

  FT _minDistSquared; // minimal distance of two points to each other
  FT _minRadiusSquared; // minimal distance of a point from center of the sphere
  FT _maxRadiusSquared; // maximal distance of a point from center of the sphere

  static constexpr bool _has_exact_rep =
    is_same_or_derived<Field_with_sqrt_tag,
                       typename Algebraic_structure_traits<FT>::Algebraic_category>::value &&
    !std::is_floating_point<FT>::value;

public:
  Delaunay_triangulation_on_sphere_traits_2(const Point_3& center = CGAL::ORIGIN,
                                            const FT radius = 1,
                                            const LK& k = LK(),
                                            const SK& sk = SK())
    : LK(k), SK(sk), _center(center), _radius(radius)
  {
    initialize_bounds();
  }

  friend void swap(Self& l, Self& r)
  {
    using std::swap;
    swap(static_cast<LK&>(l), static_cast<LK&>(r));
    swap(static_cast<SK&>(l), static_cast<SK&>(r));
    swap(l._center, r._center);
    swap(l._radius, r._radius);
    l.initialize_bounds();
    r.initialize_bounds();
  }

public:
  const LK& lk() const { return static_cast<const LK&>(*this); }
  const SK& sk() const { return static_cast<const SK&>(*this); }

  void set_center(const Point_3& center) { _center = center; }
  const Point_3& center() const { return _center; }
  void set_radius(const FT radius) { _radius = radius; initialize_bounds(); }
  FT radius() const { return _radius; }

private:
  void initialize_bounds()
  {
    Triangulations_on_sphere_2::internal::ToS2_precision_bound<FT, _has_exact_rep> bound_getter(_radius);

    _minDistSquared = bound_getter.get_squared_min_dist();
    _minRadiusSquared = bound_getter.get_squared_min_radius();
    _maxRadiusSquared = bound_getter.get_squared_max_radius();
  }

public:
  Collinear_are_strictly_ordered_on_great_circle_2
  collinear_are_strictly_ordered_on_great_circle_2_object() const
  { return Collinear_are_strictly_ordered_on_great_circle_2(center(), lk()); }

  Compare_on_sphere_2
  compare_on_sphere_2_object() const
  { return lk().compare_xyz_3_object(); }

  Equal_on_sphere_2
  equal_on_sphere_2_object() const
  { return Equal_on_sphere_2(center(), lk()); }

  Orientation_on_sphere_2
  orientation_on_sphere_2_object() const
  { return Orientation_on_sphere_2(center(), lk()); }

  Side_of_oriented_circle_on_sphere_2
  side_of_oriented_circle_on_sphere_2_object() const
  { return lk().orientation_3_object(); }

public:
  Compute_x_3
  compute_x_3_object() const
  { return lk().compute_x_3_object(); }

  Compute_y_3
  compute_y_3_object() const
  { return lk().compute_y_3_object(); }

  Compute_z_3
  compute_z_3_object() const
  { return lk().compute_z_3_object(); }

  Construct_point_on_sphere_2
  construct_point_on_sphere_2_object() const
  { return lk().construct_point_on_sphere_2_object(); }

  Construct_point_3
  construct_point_3_object() const
  { return lk().construct_point_3_object(); }

  Construct_segment_3
  construct_segment_3_object() const
  { return lk().construct_segment_3_object(); }

  Construct_triangle_3
  construct_triangle_3_object() const
  { return lk().construct_triangle_3_object(); }

  Construct_arc_on_sphere_2
  construct_arc_on_sphere_2_object() const
  { return Construct_arc_on_sphere_2(center(), radius(), lk(), sk()); }

  Construct_circumcenter_on_sphere_2
  construct_circumcenter_on_sphere_2_object() const
  { return Construct_circumcenter_on_sphere_2(center(), radius(), lk()); }

  Construct_circumcenter_3
  construct_circumcenter_3_object() const
  { return lk().construct_circumcenter_3_object();  }

public:
  bool is_on_sphere(const Point_on_sphere_2& p) const
  {
    const FT sq_dist = lk().compute_squared_distance_3_object()(p, _center);

    if(_has_exact_rep)
      return (sq_dist == _minRadiusSquared);
    else
      return (_minRadiusSquared <= sq_dist && sq_dist < _maxRadiusSquared);
  }

  bool are_points_too_close(const Point_on_sphere_2& p, const Point_on_sphere_2& q) const
  {
    if(_has_exact_rep)
      return false;
    else
      return (CGAL::squared_distance(p, q) <= _minDistSquared);
  }
};

} // namespace CGAL

#endif // CGAL_DELAUNAY_TRIANGULATION_SPHERE_TRAITS_2_H
