// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Baruch Zukerman    <baruchzu@post.tau.ac.il>
//             Efi Fogel          <efifogel@gmail.com>

#ifndef CGAL_ENV_PLANE_TRAITS_3_H
#define CGAL_ENV_PLANE_TRAITS_3_H

#include <CGAL/license/Envelope_3.h>


#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Envelope_3/Envelope_base.h>
#include <CGAL/Envelope_3/Env_plane_traits_3_functions.h>

namespace CGAL {

template <typename Kernel_,
          typename ArrLinearTraits = Arr_linear_traits_2<Kernel_>>
class Env_plane_traits_3 : public ArrLinearTraits {
public:
  using Kernel = Kernel_;
  using FT = typename Kernel::FT;
  using Base = Arr_linear_traits_2<Kernel>;
  using Self = Env_plane_traits_3<Kernel>;

  using Multiplicity = typename Base::Multiplicity;
  using Point_2 = typename Base::Point_2;
  using Curve_2 = typename Base::Curve_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
  using Plane_3 = typename Kernel::Plane_3;
  using Vector_2 = typename Kernel::Vector_2;
  using Vector_3 = typename Kernel::Vector_3;
  using Segment_2 = typename Kernel::Segment_2;
  using Ray_2 = typename Kernel::Ray_2;
  using Line_2 = typename Kernel::Line_2;
  using Line_3 = typename Kernel::Line_3;
  using Intersection_curve = std::pair<Curve_2, Multiplicity>;

  using Left_side_category = typename Base::Left_side_category;
  using Bottom_side_category = typename Base::Bottom_side_category;
  using Top_side_category = typename Base::Top_side_category;
  using Right_side_category = typename Base::Right_side_category;

  //!
  class Is_vertical_3 {
  public:
    bool operator()(const Plane_3& h) const { return CGAL::is_zero(h.c()); }
  };

  //! Obtain an Is_vertical_3 functor object.
  Is_vertical_3 is_vertical_3_object() const { return Is_vertical_3(); }

  //
  class _Env_plane {
  protected:
    Plane_3 m_plane;
    Line_2 m_line;
    bool m_is_all_plane; // true -> all plane, false -> halfplane
    bool m_is_vert;

  public:
    _Env_plane() {}

    _Env_plane(const Plane_3& h) : m_plane(h), m_is_all_plane(true) {
      Self s;
      m_is_vert = s.is_vertical_3_object()(h);
    }

    _Env_plane(const Plane_3& h, const Line_2& l) :
      m_plane(h), m_line(l), m_is_all_plane(false), m_is_vert(false) {
      CGAL_precondition_code(Self s);
      CGAL_precondition(! s.is_vertical_3_object()(h));
    }

    bool is_vertical() const { return m_is_vert; }

    const Plane_3& plane() const { return m_plane; }

    operator Plane_3() const { return m_plane; }

    const Line_2& line() const {
      CGAL_assertion(! m_is_all_plane);
      return m_line;
    }

    bool is_all_plane() const { return m_is_all_plane; }
  };

  using Xy_monotone_surface_3 = _Env_plane;
  using Surface_3 = _Env_plane;

  /*! Subdivide a given surface into \f$xy\f$-monotone parts.
   */
  class Make_xy_monotone_3 {
  public:
    template <typename OutputIterator>
    OutputIterator operator()(const Surface_3& s, bool /* is_lower */,
                              OutputIterator o) const {
      *o++ = s;
      return o;
    }
  };

  //! Obtain a Make_xy_monotone_3 functor object.
  Make_xy_monotone_3 make_xy_monotone_3_object() const
  { return Make_xy_monotone_3(); }

  /*! Determine the relative \f$z\f$-order of two given \f$xy\f$-monotone
   * surfaces at the \f$xy\f$-coordinates of a given point or \f$x\f$-monotone
   * curve.
   */
  class Compare_z_at_xy_3 {
  protected:
    using Traits_3 = Env_plane_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Compare_z_at_xy_3(const Traits_3& traits) : m_traits(traits) {}

    friend class Env_plane_traits_3<Kernel>;

  public:
    //
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      const Plane_3& plane1 = h1.plane();
      const Plane_3& plane2 = h2.plane();
      Sign sign_of_c1c2 = CGAL::sign(plane1.c() * plane2.c());
      Sign sign_of_expr =
        CGAL::sign((p.x()*plane1.a() + p.y()*plane1.b() +
                    plane1.d())*plane2.c() -
                   (p.x()*plane2.a() + p.y()*plane2.b() +
                    plane2.d())*plane1.c());
      int i = -1 * static_cast<int>(sign_of_c1c2) *
        static_cast<int>(sign_of_expr);
      return static_cast<Comparison_result>(i);
    }

    //
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      const Kernel& kernel = m_traits;
      Point_2 p = (cv.is_segment()) ?
        kernel.construct_midpoint_2_object()(cv.left(), cv.right()) :
        ((cv.is_ray()) ?
         kernel.construct_point_on_2_object()(cv.ray(), 1) :
         kernel.construct_point_on_2_object()(cv.line(), 1));

      return this->operator()(p, h1, h2);
    }

    //
    Comparison_result operator()(const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      CGAL_assertion(h1.is_all_plane() && h2.is_all_plane());

      const Plane_3& p1 = h1.plane();
      const Plane_3& p2 = h2.plane();
      const FT& res = p2.d()*p1.c() - p1.d()*p2.c();
      int i = static_cast<int>(CGAL::sign(p1.c()*p2.c())) *
        static_cast<int>(CGAL::sign (res));
      return static_cast<Comparison_result>(i);
    }
  };

  //! Obtain a Compare_z_at_xy_3 functor object.
  Compare_z_at_xy_3 compare_z_at_xy_3_object() const
  { return Compare_z_at_xy_3(*this); }

  /*! Determine the relative \f$z\f$-order of the two given \f$xy\f$-monotone
   * surfaces immediately above their projected intersection curve (a planar
   * point \f$p\f$ is above an \f$x\f$-monotone curve \f$c\f$ if it is in the
   * \f$x\f$-range of \f$c\f$, and lies to its left when the curve is traversed
   * from its \f$xy\f$-lexicographically smaller endpoint to its larger
   * endpoint).
   */
  class Compare_z_at_xy_above_3 {
  protected:
    using Traits_3 = Env_plane_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Compare_z_at_xy_above_3(const Traits_3& traits) : m_traits(traits) {}

    friend class Env_plane_traits_3<Kernel>;

  public:
    /*! Determine the relative \f$z\f$-order of the two given \f$xy\f$-monotone
     * surfaces immediately above their projected intersection curve, which is
     * also given.
     *
     * \param cv the intersection curve.
     * \param h1 the first surface.
     * \param h2 the second surface.
     * \pre `h1` and `h2` are defined "above" `cv`, and their relative
     * \f$z\f$-order is the same for some small enough neighborhood of points
     * above `cv`.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      const Plane_3& plane1 = h1.plane();
      const Plane_3& plane2 = h2.plane();

      const FT& a1 = plane1.a();
      const FT& b1 = plane1.b();
      const FT& c1 = plane1.c();

      const FT& a2 = plane2.a();
      const FT& b2 = plane2.b();
      const FT& c2 = plane2.c();

      // our line is a3*x + b3*y + c3 = 0
      // it is assumed that the planes intersect over this line
      const Line_2& line = cv.supp_line();
      const FT& a3 = line.a();
      const FT& b3 = line.b();
      // const FT& c3 = line.c();       // unused

      // if the line was parallel to the y-axis (i.e x = const),
      // then it was enough to compare dz/dx of both planes
      // for general line, we change coordinates to (v, w), preserving
      // orientation, so the line is the w-axis in the new coordinates
      // (i.e v = const).
      //
      // ( v )  =  A ( x )    where A = (  a3  b3 )
      //   w           y                  -b3  a3
      //
      // so v =  a3*x + b3*y
      //    w = -b3*x + a3*y
      // preserving orientation since detA = a3^2 +b3^2 > 0
      //
      // We compute the planes equations in the new coordinates
      // and compare dz/dv
      //
      // ( x )  =  A^(-1) ( v )    where A^(-1) = ( a3  -b3 ) * detA^(-1)
      //   y                w                       b3   a3
      // so x = (a3*v - b3*w)*(1/detA)
      //    y = (b3*v + a3*w)*(1/detA)
      // plane1 ==> (a1a3 + b1b3)v + (b1a3 - a1b3)w + (c1z + d1)*detA = 0
      // plane2 ==> (a2a3 + b2b3)v + (b2a3 - a2b3)w + (c2z + d2)*detA = 0
      //
      // dz/dv(1) = (-a1a3 - b1b3) / c1*detA
      // dz/dv(2) = (-a2a3 - b2b3) / c2*detA
      // since detA>0 we can omit it.
      //
      Sign s1 = CGAL_NTS sign((a2*a3+b2*b3)/c2-(a1*a3+b1*b3)/c1);

      // We only need to make sure that w is in the correct direction
      // (going from down to up)
      // the original segment endpoints p1=(x1,y1) and p2=(x2,y2)
      // are transformed to (v1,w1) and (v2,w2), so we need that w2 > w1
      // (otherwise the result should be multiplied by -1)

      const Kernel& kernel = m_traits;
      Point_2 p1 = kernel.construct_point_on_2_object()(line, 0);
      Point_2 p2 = kernel.construct_point_on_2_object()(line, 1);

      if (kernel.compare_xy_2_object()(p1, p2) == LARGER) std::swap(p1, p2);

      CGAL_assertion(kernel.compare_xy_2_object()(p1, p2) == SMALLER);

      const FT& x1 = p1.x();
      const FT& y1 = p1.y();
      const FT& x2 = p2.x();
      const FT& y2 = p2.y();

      Sign s2 = CGAL_NTS sign(-b3*x1+a3*y1-(-b3*x2+a3*y2));
      return s1 * s2;
    }
  };

  //! Obtain a Compare_z_at_xy_above_3 functor object.
  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object() const
  { return Compare_z_at_xy_above_3(*this); }

  /*! Determine the relative \f$z\f$-order of the two given \f$xy\f$-monotone
   * surfaces immediately below their projected intersection curve (a planar
   * point \f$p\f$ is below an \f$x\f$-monotone curve \f$c\f$ if it is in the
   * \f$x\f$-range of \f$c\f$, and lies to its left when the curve is traversed
   * from its \f$xy\f$-lexicographically smaller endpoint to its larger
   * endpoint).
   */
  class Compare_z_at_xy_below_3 {
  protected:
    using Traits_3 = Env_plane_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Compare_z_at_xy_below_3(const Traits_3& traits) : m_traits(traits) {}

    friend class Env_plane_traits_3<Kernel>;

  public:
    /*! Determine the relative \f$z\f$-order of the two given \f$xy\f$-monotone
     * surfaces immediately below their projected intersection curve, which is
     * also given.
     *
     * \param cv the intersection curve.
     * \param h1 the first surface.
     * \param h2 the second surface.
     * \pre `h1` and `h2` are defined "above" `cv`, and their relative
     * \f$z\f$-order is the same for some small enough neighborhood of points
     * below `cv`.
     */
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      auto cmp_above = m_traits.compare_z_at_xy_above_3_object();
      return CGAL::opposite(cmp_above(cv, h1, h2));
    }
  };

  //! Obtain a Compare_z_at_xy_below_3 functor object.
  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object() const
  { return Compare_z_at_xy_below_3(*this); }

  /*! Compute all planar \f$x\f$-monotone curves and possibly isolated planar
   * points that form the projection of the boundary of the given
   * \f$xy\f$-monotone surface s onto the \f$xy\f$-plane.
   */
  class Construct_projected_boundary_2 {
  protected:
    using Traits_3 = Env_plane_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Construct_projected_boundary_2(const Traits_3& traits) : m_traits(traits) {}

    friend class Env_plane_traits_3<Kernel>;

  public:
    template <typename OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s,
                              OutputIterator o) const {
      if (s.is_all_plane()) {
        if (! s.is_vertical()) return o;

        const Plane_3& h = s.plane();
        Line_2 proj_line(h.a(), h.b(), h.d());
        *o++ = std::make_pair(X_monotone_curve_2(proj_line),
                              ON_ORIENTED_BOUNDARY);
        return o;
      }

      // s is half-plane
      const Kernel& kernel = m_traits;
      const Point_2& p1 = kernel.construct_point_on_2_object()(s.line(), 0);
      const Point_2& p2 = kernel.construct_point_on_2_object()(s.line(), 1);
      Comparison_result res = kernel.compare_xy_2_object()(p1, p2);

      Oriented_side side =
        (res == SMALLER) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
      *o++ = std::make_pair(X_monotone_curve_2(s.line()), side);
      return o;
    }
  };

  //
  Construct_projected_boundary_2
  construct_projected_boundary_2_object() const
  { return Construct_projected_boundary_2(*this); }

  /*! compute the projection of the intersections of the \f$xy\f$-monotone
   * surfaces onto the \f$xy\f$-plane,
   */
  class Construct_projected_intersections_2 {
  protected:
    using Traits_3 = Env_plane_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Construct_projected_intersections_2(const Traits_3& traits) :
      m_traits(traits)
    {}

    friend class Env_plane_traits_3<Kernel>;

  public:
    template <typename OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s1,
                              const Xy_monotone_surface_3& s2,
                              OutputIterator o) const {
      const Kernel& kernel = m_traits;

      const Plane_3& h1 = s1.plane();
      const Plane_3& h2 = s2.plane();

      if (s1.is_vertical() && s2.is_vertical()) {
        Line_2 l1(h1.a(), h1.b(), h1.d());
        Line_2 l2(h2.a(), h2.b(), h2.d());
        auto obj = kernel.intersect_2_object()(l1, l2);

        if (const auto* p = std::get_if<Point_2>(&(*obj))) *o++ = *p;

        // otherwise, the vertical planes are parallel or overlap, so we return
        // nothing.
        return o;
      }

      if (s1.is_all_plane() && s2.is_all_plane()) {
        auto obj = kernel.intersect_3_object()(h1, h2);
        CGAL_assertion(obj != std::nullopt);
        if (const auto* l = std::get_if<Line_3>(&(*obj)))
          *o++ = Intersection_curve(project_xy(*l, kernel), 1);
        return o;
      }

      if (s1.is_all_plane() && ! s2.is_all_plane()) {
        auto obj = plane_half_plane_proj_intersection(h1, h2, s2.line(), kernel);
        if (obj == std::nullopt) return o;
        if (const auto* line = std::get_if<Line_2>(&(*obj))) {
          *o++ = Intersection_curve(*line, 1);
          return o;
        }
        if (const auto* ray = std::get_if<Ray_2>(&(*obj))) {
          *o++ = Intersection_curve(*ray, 1);
          return o;
        }
        return o;
      }
      if (! s2.is_all_plane() && s2.is_all_plane()) {
        auto obj = plane_half_plane_proj_intersection(h2, h1, s1.line(), kernel);
        if (obj == std::nullopt) return o;
        if (const auto* line = std::get_if<Line_2>(&(*obj))) {
          *o++ = Intersection_curve(*line, 1);
          return o;
        }
        if (const auto* ray = std::get_if<Ray_2>(&(*obj))) {
          *o++ = Intersection_curve(*ray, 1);
          return o;
        }
        return o;
      }

      CGAL_assertion(! s2.is_all_plane() && ! s2.is_all_plane());
      auto obj = half_plane_half_plane_proj_intersection(h1, s1.line(),
                                                         h2, s2.line(), kernel);

      if (obj == std::nullopt) return o;

      if (const auto* line = std::get_if<Line_2>(&(*obj))) {
        *o++ = Intersection_curve(*line, 1);
        return o;
      }
      if (const auto* ray = std::get_if<Ray_2>(&(*obj))) {
        *o++ = Intersection_curve(*ray, 1);
        return o;
      }
      if (const auto* seg = std::get_if<Segment_2>(&(*obj))) {
        *o++ = Intersection_curve(*seg, 1);
        return o;
      }
      if (const auto* p = std::get_if<Point_2>(&(*obj))) {
        *o++ = *p;
        return o;
      }
      return o;
    }
  };

  //! Obtain a Construct_projected_intersections_2 functor object.
  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const
  { return Construct_projected_intersections_2(*this); }
};

} //namespace CGAL

#endif
