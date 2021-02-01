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

#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Has_conversion.h>
#include <CGAL/Spherical_kernel_3.h>

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
    if(collinear(_center, p, r) || // @todo use LK
       collinear(_center, q, r) ||
       orientation(_center, p, q, r) != COLLINEAR)
    {
      return false;
    }

    if(collinear(_center, p, q))
      return true;

    const Orientation op = _lk.coplanar_orientation_3_object()(_center, p, q, r);
    const Orientation oq = _lk.coplanar_orientation_3_object()(_center, q, p, r); // @todo add an early exit
    CGAL_assertion(op != COLLINEAR && oq != COLLINEAR);

    return (op == POSITIVE) && (oq == POSITIVE);
  }

protected:
  const Point_3& _center;
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
    // the orientation of the plane (and thus of the circle) does not really matter
    // since the construction of the circular arc uses the positive normal of the circle's
    // supporting plane...
    typedef typename CGAL::internal::Converter_selector<LK, SK>::type Converter;
    Converter cv;

    typename SK::Point_3 sc = cv(_center); // @tmp cache that, I guess
    typename SK::Point_3 sp = cv(p);
    typename SK::Point_3 sq = cv(q);

    typename SK::Plane_3 pl = _sk.construct_plane_3_object()(sp, sc, sq);
    typename SK::Circle_3 c = _sk.construct_circle_3_object()(sc, square(_radius), pl);

    typename SK::Circular_arc_point_3 cp = _sk.construct_circular_arc_point_3_object()(sp);
    typename SK::Circular_arc_point_3 cq = _sk.construct_circular_arc_point_3_object()(sq);

    // @fixme ensure in arc_dual and arc_segment that 'cp' and 'cq' are in the correct order
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
//  : public LK_, // @tmp disabled to see what is really needed
//    public SK_
{
  typedef LK_                                                        LK;
  typedef SK_                                                        SK;
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

  // For the hilbert sort, not actually needed when the traits derives from LK @todo
  typedef typename LK::Compute_x_3                                   Compute_x_3;
  typedef typename LK::Compute_y_3                                   Compute_y_3;
  typedef typename LK::Compute_z_3                                   Compute_z_3;
  Compute_x_3 compute_x_3_object() const { return Compute_x_3(); }
  Compute_y_3 compute_y_3_object() const { return Compute_y_3(); }
  Compute_z_3 compute_z_3_object() const { return Compute_z_3(); }

  typedef internal::Construct_arc_on_sphere_2<LK, SK>                Construct_arc_on_sphere_2;

  // @todo needs to offer dual_on_sphere that projects on the sphere
  typedef void Construct_circumcenter_on_sphere_2;

  typedef typename LK::Construct_circumcenter_3                      Construct_circumcenter_3;

  Delaunay_triangulation_on_sphere_traits_2(const Point_3& center = CGAL::ORIGIN,
                                            const FT radius = 1,
                                            const LK& k = LK(),
                                            const SK& sk = SK())
    : /*LK(k), SK(sk),*/ _center(center), _radius(radius)
  {
    initialize_bounds();
  }

private:
  void initialize_bounds()
  {
    if(_has_exact_rep)
    {
      _minDistSquared = 0;
      _minRadiusSquared = _maxRadiusSquared = square(_radius);
    }
    else
    {
      // @fixme
      // - check the correctness of the implementation (not all FT have double precision)
      // - if LK can represent algebraic coordinates, there is no need for that
      const FT minDist = _radius * std::pow(2, -23);
      const FT minRadius = _radius * (1 - std::pow(2, -50));
      const FT maxRadius = _radius * (1 + std::pow(2, -50));
      _minDistSquared = CGAL::square(minDist);
      _minRadiusSquared = CGAL::square(minRadius);
      _maxRadiusSquared = CGAL::square(maxRadius);
    }
  }

public:
  Collinear_are_strictly_ordered_on_great_circle_2
  collinear_are_strictly_ordered_on_great_circle_2_object() const
  { return Collinear_are_strictly_ordered_on_great_circle_2(_center, LK()); }

  Compare_on_sphere_2
  compare_on_sphere_2_object() const
  { return typename LK::Compare_xyz_3(); } // @tmp static_cast<const Base&>(*this)

  Equal_on_sphere_2
  equal_on_sphere_2_object() const
  { return Equal_on_sphere_2(_center, LK()); } // @tmp

  Orientation_on_sphere_2
  orientation_on_sphere_2_object() const
  { return Orientation_on_sphere_2(_center, LK()); }

  Side_of_oriented_circle_on_sphere_2
  side_of_oriented_circle_on_sphere_2_object() const
  { return typename LK::Orientation_3(); } // @tmp

public:
  Construct_point_on_sphere_2
  construct_point_on_sphere_2_object() const
  { return Construct_point_on_sphere_2(); } // @tmp

  Construct_point_3
  construct_point_3_object() const
  { return typename LK::Construct_point_3(); } // @tmp

  Construct_segment_3
  construct_segment_3_object() const
  { return typename LK::construct_segment_3_object(); } // @tmp

  Construct_triangle_3
  construct_triangle_3_object() const
  { return typename LK::Construct_triangle_3(); } // @tmp

  Construct_arc_on_sphere_2
  construct_arc_on_sphere_2_object() const
  { return Construct_arc_on_sphere_2(_center, _radius, LK(), SK()); } // @tmp

  Construct_circumcenter_on_sphere_2
  construct_circumcenter_on_sphere_2_object() const
  { return typename LK::Construct_circumcenter_3(); } // @tmp

  Construct_circumcenter_3
  construct_circumcenter_3_object() const
  { return typename LK::Construct_circumcenter_3();  } // @tmp

public:
  void set_center(const Point_3& center) { _center = center; }
  const Point_3& center() const { return _center; }
  void set_radius(const FT radius) { _radius = radius; initialize_bounds(); }
  FT radius() const { return _radius; }

public:
  bool is_on_sphere(const Point_on_sphere_2& p) const
  {
    const FT sq_dist = CGAL::squared_distance(p, _center); // @todo use LK

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
};

} // namespace CGAL

#endif // CGAL_DELAUNAY_TRIANGULATION_SPHERE_TRAITS_2_H
