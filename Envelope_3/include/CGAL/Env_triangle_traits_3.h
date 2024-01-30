// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//             Baruch Zukerman        <baruchzu@post.tau.ac.il>
//             Efi Fogel              <efifogel@gmail.com>

/*! \file CGAL/Envelope_triangles_traits_3.h
 * \brief Model for CGAL's EnvelopeTraits_3 concept.
 * \endlink
 */

#ifndef CGAL_ENV_TRIANGLE_TRAITS_3_H
#define CGAL_ENV_TRIANGLE_TRAITS_3_H

#include <CGAL/license/Envelope_3.h>


#include <vector>

#include <CGAL/enum.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Envelope_3/Envelope_base.h>

namespace CGAL {

template <typename Kernel_> class Env_triangle_3;

// this traits class supports both triagles and segments in 3d
template <typename Kernel_,
          typename ArrSegmentTraits = Arr_segment_traits_2<Kernel_>>
class Env_triangle_traits_3 : public ArrSegmentTraits {
public:
  using Traits_2 = Arr_segment_traits_2<Kernel_>;
  using Point_2 = typename Traits_2::Point_2;
  using X_monotone_curve_2 = typename Traits_2::X_monotone_curve_2;
  using Multiplicity = typename Traits_2::Multiplicity;

  using Kernel = Kernel_;
  using Self = Env_triangle_traits_3<Kernel>;

  using Point_3 = typename Kernel::Point_3;

  /*! \class Representation of a 3d triangle with cached data.
   */
  class _Triangle_cached_3 {
  public:
    using Plane_3 = typename Kernel::Plane_3;
    using Triangle_3 = typename Kernel::Triangle_3;
    using Point_3 = typename Kernel::Point_3;
    using Segment_3 = typename Kernel::Segment_3;

  protected:
    Plane_3 pl;          // The plane that supports the triangle.
    Point_3 vertices[3]; // The vertices of the triangle.
    bool is_vert;     // Is this a vertical triangle (or a segment).
    bool is_seg;  // Is this a segment.

  public:
    /*! Default constructor.
     */
    _Triangle_cached_3() : is_vert(false), is_seg(false) {}

    /*! Constructor from a non-degenerate triangle.
     * \param tri The triangle.
     * \pre The triangle is not degenerate.
     */
    _Triangle_cached_3(const Triangle_3& tri) {
      Kernel kernel;
      CGAL_assertion(! kernel.is_degenerate_3_object()(tri));

      auto construct_vertex = kernel.construct_vertex_3_object();

      vertices[0] = construct_vertex(tri, 0);
      vertices[1] = construct_vertex(tri, 1);
      vertices[2] = construct_vertex(tri, 2);

      pl = kernel.construct_plane_3_object()(vertices[0],
                                             vertices[1], vertices[2]);
      Self self;
      is_vert = kernel.collinear_2_object()(self.project(vertices[0]),
                                            self.project(vertices[1]),
                                            self.project(vertices[2]));
      is_seg = false;
    }

    /*! Construct a triangle from three non-collinear end-points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \param p3 The third point.
     * \pre The 3 endpoints are not the collinear.
     */
    _Triangle_cached_3(const Point_3& p1, const Point_3& p2,
                       const Point_3& p3) {
      Kernel kernel;
      CGAL_assertion(! kernel.collinear_3_object()(p1, p2, p3));

      vertices[0] = p1;
      vertices[1] = p2;
      vertices[2] = p3;

      pl = kernel.construct_plane_3_object()(vertices[0], vertices[1],
                                             vertices[2]);
      Self self;
      is_vert = kernel.collinear_2_object()(self.project(vertices[0]),
                                            self.project(vertices[1]),
                                            self.project(vertices[2]));
      is_seg = false;
    }

    /*! Construct a triangle from 3 end-points on a supporting plane.
     * \param supp_plane The supporting plane.
     * \param p1 The first point.
     * \param p2 The second point.
     * \param p3 The third point.
     * \pre The 3 endpoints are not the collinear and all lie on the given
     *      plane.
     */
    _Triangle_cached_3(const Plane_3& supp_plane, const Point_3& p1,
                       const Point_3& p2, const Point_3& p3) :
      pl(supp_plane) {
      Kernel kernel;

      CGAL_precondition(kernel.has_on_3_object() (pl, p1) &&
                        kernel.has_on_3_object() (pl, p2) &&
                        kernel.has_on_3_object() (pl, p3));
      CGAL_precondition(!kernel.collinear_3_object()(p1, p2, p3));

      vertices[0] = p1;
      vertices[1] = p2;
      vertices[2] = p3;

      Self self;
      is_vert = kernel.collinear_2_object()(self.project(vertices[0]),
                                            self.project(vertices[1]),
                                            self.project(vertices[2]));
      is_seg = false;
    }

    /*! Constructor from a segment.
     * \param seg The segment.
     * \pre The segment is not degenerate.
     */
    _Triangle_cached_3(const Segment_3& seg) {
      Kernel kernel;
      CGAL_assertion(! kernel.is_degenerate_3_object()(seg));

      typename Kernel::Construct_vertex_3
        construct_vertex = kernel.construct_vertex_3_object();

      vertices[0] = construct_vertex(seg, 0);
      vertices[1] = construct_vertex(seg, 1);
      vertices[2] = vertices[1];

      is_vert = true;
      is_seg = true;

      // construct a vertical plane through the segment
      Point_3 tmp(vertices[0].x(), vertices[0].y(), vertices[0].z()-1);
      pl = kernel.construct_plane_3_object()(vertices[0], vertices[1], tmp);
    }

    /*! Constructor from two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \param seg The segment.
     * \pre The segment between the points is not degenerate.
     */
    _Triangle_cached_3(const Point_3& p1, const Point_3& p2) {
      Kernel kernel;
      CGAL_assertion(! kernel.equal_3_object()(p1, p2));

      vertices[0] = p1;
      vertices[1] = p2;
      vertices[2] = p2;

      is_vert = true;
      is_seg = true;

      // construct a vertical plane through the segment
      Point_3 tmp(vertices[0].x(), vertices[0].y(), vertices[0].z()-1);
      pl = kernel.construct_plane_3_object()(vertices[0], vertices[1], tmp);
    }

    /*! Assignment operator.
     * \param tri the source triangle to copy from
     */
    const _Triangle_cached_3& operator=(const Triangle_3& tri) {
      Kernel kernel;
      CGAL_assertion(! kernel.is_degenerate_3_object()(tri));

      auto construct_vertex = kernel.construct_vertex_3_object();

      vertices[0] = construct_vertex(tri, 0);
      vertices[1] = construct_vertex(tri, 1);
      vertices[2] = construct_vertex(tri, 2);

      pl = kernel.construct_plane_3_object()(vertices[0],
                                             vertices[1], vertices[2]);
      Self self;
      is_vert = kernel.collinear_2_object()(self.project(vertices[0]),
                                            self.project(vertices[1]),
                                            self.project(vertices[2]));
      is_seg = false;

      return (*this);
    }

    /*! Get the ith endpoint.
     */
    const Point_3& vertex(unsigned int i) const { return vertices[i%3]; }

    /*! Get the supporting plane.
     */
    const Plane_3& plane() const { return (pl); }

    /*! Check whether the triangle is vertical.
     */
    bool is_vertical() const { return (is_vert); }

    /*! Check whether the surface is a segment.
     */
    bool is_segment() const { return (is_seg); }

    /*! Check whether the surface is xy-monotone (false, if it is a vertical
     * triangle)
     */
    bool is_xy_monotone() const { return (!is_vertical() || is_segment()); }
  };

public:
  // types for EnvelopeTraits_3 concept
  //! type of xy-monotone surfaces
  using Xy_monotone_surface_3 = Env_triangle_3<Kernel>;
  //! type of surfaces
  using Surface_3 = Xy_monotone_surface_3;

  // we have a collision between the Kernel's Intersect_2 and the one
  // from the segment traits
  using Intersect_2 = typename Traits_2::Intersect_2;

protected:
  using FT = typename Kernel::FT;
  using Triangle_2 = typename Kernel::Triangle_2;
  using Segment_2 = typename Kernel::Segment_2;

  using Segment_3 = typename Kernel::Segment_3;
  using Triangle_3 = typename Kernel::Triangle_3;
  using Plane_3 = typename Kernel::Plane_3;

  using Assign_2 = typename Kernel::Assign_2;
  using Construct_vertex_2 = typename Kernel::Construct_vertex_2;

  using Assign_3 = typename Kernel::Assign_3;
  using Intersect_3 = typename Kernel::Intersect_3;
  using Construct_vertex_3 = typename Kernel::Construct_vertex_3;


  using Line_2 = typename Kernel::Line_2;
  using Direction_2 = typename Kernel::Direction_2;

  using Line_3 = typename Kernel::Line_3;
  using Direction_3 = typename Kernel::Direction_3;

  using Intersection_curve = std::pair<X_monotone_curve_2, Multiplicity>;

public:
  /***************************************************************************/
  // EnvelopeTraits_3 functors
  /***************************************************************************/

  /*! Subdivide a given surface into \f$xy\f$-monotone parts.
   */
  class Make_xy_monotone_3 {
  protected:
    using Traits_3 = Env_triangle_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Make_xy_monotone_3(const Traits_3& traits) : m_traits(traits) {}

    friend class Env_triangle_traits_3<Kernel>;

  public:
    // create xy-monotone surfaces from a general surface
    // return a past-the-end iterator
    template <typename OutputIterator>
    OutputIterator operator()(const Surface_3& s, bool is_lower,
                              OutputIterator o) const {
      m_is_lower = is_lower;

      // a non-vertical triangle is already xy-monotone
      if (! s.is_vertical()) *o++ = s;
      else {
        // split a vertical triangle into one or two segments
        const Point_3& a1 = s.vertex(0);
        const Point_3& a2 = s.vertex(1);
        const Point_3& a3 = s.vertex(2);
        Point_2 b1 = m_traits.project(a1);
        Point_2 b2 = m_traits.project(a2);
        Point_2 b3 = m_traits.project(a3);
        const Kernel& k = m_traits;
        if (k.collinear_are_ordered_along_line_2_object()(b1, b2, b3)) {
          if (k.equal_2_object()(b1, b2))
            // only one segment in the output - the vertical does not count
            *o++ = Xy_monotone_surface_3(find_envelope_point(a1, a2), a3);
          else if (k.equal_2_object()(b2, b3))
            *o++ = Xy_monotone_surface_3(a1, find_envelope_point(a2, a3));
          else
            // check whether two or one segments appear on the envelope
            return find_envelope_segments(a1, a2, a3, s.plane(), o);
        }
        else if (k.collinear_are_ordered_along_line_2_object()(b1, b3, b2)) {
          if (k.equal_2_object()(b1, b3))
            // only one segment in the output
            *o++ = Xy_monotone_surface_3(find_envelope_point(a1, a3), a2);
          else
            // check whether two or one segments appear on the envelope
            return find_envelope_segments(a1, a3, a2, s.plane(), o);
        }
        else {
          // check whether two or one segments appear on the envelope
          return find_envelope_segments(a2, a1, a3, s.plane(), o);
        }
      }
      return o;
    }

  protected:
    // find the envelope point among the two points with same xy coordinates
    const Point_3& find_envelope_point(const Point_3& p1, const Point_3& p2)
      const {
      CGAL_precondition(p1.x() == p2.x() && p1.y() == p2.y());
      const Kernel& k = m_traits;
      Comparison_result cr = k.compare_z_3_object()(p1, p2);
      CGAL_assertion(cr != EQUAL);
      if ((m_is_lower && cr == SMALLER) || (! m_is_lower && cr == LARGER))
        return p1;
      else return p2;
    }

    // get the three triangle coordinates (ordered along 2d-line) and find
    // which segment(s) is(are) the envelope of this triangle
    // "plane" is the vertical plane on which the triangle lies
    template <typename OutputIterator>
    OutputIterator find_envelope_segments(const Point_3& p1, const Point_3& p2,
                                          const Point_3& p3,
                                          const Plane_3& plane,
                                          OutputIterator o) const {
      // our vertical plane is a*x + b*y + d = 0
      FT a = plane.a(), b = plane.b();
      CGAL_precondition(plane.c() == 0);

      // if the plane was parallel to the yz-plane (i.e x = const),
      // then it was enough to use the y,z coordinates as in the 2-dimensional
      // case, to find whether a 2d point lies below/above/on a line
      // this test is simply computing the sign of:
      //    (1)    [(y3 - y1)(z2 - z1) - (z3 - z1)(y2 - y1)] * sign(y3 - y1)
      // abd comparing to 0, where pi = (xi, yi, zi), and p2 is compared to the
      // line formed by p1 and p3 (in the direction p1 -> p3)
      //
      // for general vertical plane, we change (x, y) coordinates to (v, w),
      // (keeping the z-coordinate as is)
      // so the plane is parallel to the wz-plane in the new coordinates
      // (i.e v = const).
      //
      // ( v )  =  A ( x )    where A = (  a  b )
      //   w           y                  -b  a
      //
      // so v =  a*x + b*y
      //    w = -b*x + a*y
      //
      // Putting the new points coordinates in equation (1) we get:
      //    (2)    (w3 - w1)(z2 - z1) - (z3 - z1)(w2 - w1) =
      //           (-b*x3 + a*y3 + b*x1 - a*y1)(z2 - z1) -
      //                                 (z3 - z1)(-b*x2 + a*y2 + b*x1 - a*y1)
      //
      FT w1 = a*p1.y() - b*p1.x();
      FT w2 = a*p2.y() - b*p2.x();
      FT w3 = a*p3.y() - b*p3.x();

      Sign s1 = CGAL::sign((w3 - w1)*(p2.z() - p1.z()) -
                           (p3.z() - p1.z())*(w2 - w1));

      // the points should not be collinear
      CGAL_assertion(s1 != 0);

      // should also take care for the original and transformed direction of
      // the segment
      Sign s2 = CGAL_NTS sign(w3 - w1);
      Sign s = CGAL_NTS sign(int(s1 * s2));

      bool use_one_segment = true;
      if ((m_is_lower && (s == NEGATIVE)) || (! m_is_lower && (s == POSITIVE)))
        use_one_segment = false;

      if (use_one_segment) *o++ = Xy_monotone_surface_3(p1, p3);
      else {
        *o++ = Xy_monotone_surface_3(p1, p2);
        *o++ = Xy_monotone_surface_3(p2, p3);
      }
      return o;
    }

    mutable bool m_is_lower;
  };

  /*! Obtain a Make_xy_monotone_3 functor object. */
  Make_xy_monotone_3
  make_xy_monotone_3_object() const { return Make_xy_monotone_3(*this); }

  /*! Compute all planar \f$x\f$-monotone curves and possibly isolated planar
   * points that form the projection of the boundary of the given
   * \f$xy\f$-monotone surface s onto the \f$xy\f$-plane.
   */
  class Construct_projected_boundary_2 {
  protected:
    using Traits_3 = Env_triangle_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Construct_projected_boundary_2(const Traits_3& traits) : m_traits(traits) {}

    friend class Env_triangle_traits_3<Kernel>;

  public:
    // insert into the OutputIterator all the (2d) curves of the boundary of
    // the vertical projection of the surface on the xy-plane
    // the OutputIterator value type is X_monotone_curve_2
    template <typename OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s,
                              OutputIterator o) const {
      // the input xy-monotone surface should be either non-vertical or
      // a segment
      CGAL_assertion(s.is_xy_monotone());

      if (! s.is_vertical()) {
        // the projection is a triangle
        const Point_3& a1 = s.vertex(0);
        const Point_3& a2 = s.vertex(1);
        const Point_3& a3 = s.vertex(2);
        Point_2 b1 = m_traits.project(a1);
        Point_2 b2 = m_traits.project(a2);
        Point_2 b3 = m_traits.project(a3);

        const Kernel& k = m_traits;

        X_monotone_curve_2 A(b1, b2);
        X_monotone_curve_2 B(b2, b3);
        X_monotone_curve_2 C(b3, b1);

        const Line_2& l1 =
          (A.is_directed_right()) ? A.line() : A.line().opposite();
        const Line_2& l2 =
          (B.is_directed_right()) ? B.line() : B.line().opposite();
        const Line_2& l3 =
          (C.is_directed_right()) ? C.line() : C.line().opposite();

        Oriented_side s1 = k.oriented_side_2_object()(l1, b3);
        Oriented_side s2 = k.oriented_side_2_object()(l2, b1);
        Oriented_side s3 = k.oriented_side_2_object()(l3, b2);

        CGAL_assertion((s1 != ON_ORIENTED_BOUNDARY) &&
                       (s2 != ON_ORIENTED_BOUNDARY) &&
                       (s3 != ON_ORIENTED_BOUNDARY));

        *o++ = std::make_pair(A, s1);
        *o++ = std::make_pair(B, s2);
        *o++ = std::make_pair(C, s3);
        return o;
      }

      // s is a segment, and so is its projection
      // s shouldn't be a z-vertical segment
      const Point_3& a1 = s.vertex(0);
      const Point_3& a2 = s.vertex(1);

      Point_2 b1 = m_traits.project(a1);
      Point_2 b2 = m_traits.project(a2);
      CGAL_assertion(b1 != b2);

      *o++ = std::make_pair(X_monotone_curve_2(b1, b2), ON_ORIENTED_BOUNDARY);
      return o;
    }
  };

  /*! Obtain a Construct_projected_boundary_curves_2 functor object. */
  Construct_projected_boundary_2
  construct_projected_boundary_2_object() const
  { return Construct_projected_boundary_2(*this); }

  /*! compute the projection of the intersections of the \f$xy\f$-monotone
   * surfaces onto the \f$xy\f$-plane,
   */
  class Construct_projected_intersections_2 {
  protected:
    using Traits_3 = Env_triangle_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Construct_projected_intersections_2(const Traits_3& traits) :
      m_traits(traits)
    {}

    friend class Env_triangle_traits_3<Kernel>;

  public:
    // insert into OutputIterator all the (2d) projections on the xy plane of
    // the intersection objects between the 2 surfaces
    // the data type of OutputIterator is Object
    template <typename OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s1,
                              const Xy_monotone_surface_3& s2,
                              OutputIterator o) const {
      CGAL_assertion(s1.is_xy_monotone() && s2.is_xy_monotone());

      if (! m_traits.do_intersect(s1, s2)) return o;

      const Kernel& k = m_traits;
      std::optional<std::variant<Point_3, Segment_3>> inter_obj =
        m_traits.intersection(s1, s2);
      if (inter_obj == std::nullopt) return o;

      if (const auto* point = std::get_if<Point_3>(&(*inter_obj))) {
        *o++ = m_traits.project(*point);
        return o;
      }

      const auto* curve = std::get_if<Segment_3>(&(*inter_obj));
      CGAL_assertion(curve != nullptr);

      Segment_2 proj_seg = m_traits.project(*curve);
      if (! k.is_degenerate_2_object() (proj_seg)) {
        Intersection_curve inter_cv (proj_seg, 1);
        *o++ = inter_cv;
        return o;
      }

      const Point_2& p = k.construct_point_on_2_object() (proj_seg, 0);
      *o++ = p;
      return o;
    }
  };

  /*! Obtain a Construct_projected_intersections_2 functor object. */
  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const
  { return Construct_projected_intersections_2(*this); }

  /*! Determine the relative \f$z\f$-order of two given \f$xy\f$-monotone
   * surfaces at the \f$xy\f$-coordinates of a given point or \f$x\f$-monotone
   * curve.
   */
  class Compare_z_at_xy_3 {
  protected:
    using Traits_3 = Env_triangle_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Compare_z_at_xy_3(const Traits_3& traits) : m_traits(traits) {}

    friend class Env_triangle_traits_3<Kernel>;

  public:
    // check which of the surfaces is closer to the envelope at the xy
    // coordinates of point
    // (i.e. lower if computing the lower envelope, or upper if computing
    // the upper envelope)
    // precondition: the surfaces are defined in point
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& surf1,
                                 const Xy_monotone_surface_3& surf2) const {
      // we compute the points on the planes, and then compare their z
      // coordinates
      const Plane_3& plane1 = surf1.plane();
      const Plane_3& plane2 = surf2.plane();

      // if the 2 triangles have the same supporting plane, and they are not
      // vertical, then they have the same z coordinate over this point
      if ((plane1 == plane2 || plane1 == plane2.opposite()) &&
          ! surf1.is_vertical())
        return EQUAL;

      // Compute the intersetion between the vertical line and the given
      // surfaces
      const Kernel& k = m_traits;
      Point_3 ip1 = m_traits.envelope_point_of_surface(p, surf1);
      Point_3 ip2 = m_traits.envelope_point_of_surface(p, surf2);
      return k.compare_z_3_object()(ip1, ip2);
    }

    // check which of the surfaces is closer to the envelope at the xy
    // coordinates of cv
    // (i.e. lower if computing the lower envelope, or upper if computing the
    // upper envelope)
    // precondition: the surfaces are defined in all points of cv,
    //               and the answer is the same for each of these points
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& surf1,
                                 const Xy_monotone_surface_3& surf2) const {
      // first try the endpoints, if cannot be sure, use the mid point
      Comparison_result res =
        m_traits.compare_z_at_xy_3_object()(cv.left(), surf1, surf2);

      if (res == EQUAL) {
        res = m_traits.compare_z_at_xy_3_object()(cv.right(), surf1, surf2);
        if (res == EQUAL) {
          Point_2 mid = m_traits.construct_middle_point(cv);
          res = m_traits.compare_z_at_xy_3_object()(mid, surf1, surf2);
        }
      }

      return res;
    }
  };

  /*! Obtain a Compare_z_at_xy_3 functor object. */
  Compare_z_at_xy_3
  compare_z_at_xy_3_object() const { return Compare_z_at_xy_3(*this); }

  /*! Determine the relative \f$z\f$-order of the two given \f$xy\f$-monotone
   * surfaces immediately above their projected intersection curve (a planar
   * point \f$p\f$ is above an \f$x\f$-monotone curve \f$c\f$ if it is in the
   * \f$x\f$-range of \f$c\f$, and lies to its left when the curve is traversed
   * from its \f$xy\f$-lexicographically smaller endpoint to its larger
   * endpoint).
   */
  class Compare_z_at_xy_above_3 {
  protected:
    using Traits_3 = Env_triangle_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Compare_z_at_xy_above_3(const Traits_3& traits) : m_traits(traits) {}

    friend class Env_triangle_traits_3<Kernel>;

  public:
    // check which of the surfaces is closer to the envelope on the points
    // above the curve cv
    // (i.e. lower if computing the lower envelope, or upper if computing the
    // upper envelope)
    // precondition: the surfaces are defined above cv (to the left of cv,
    //               if cv is directed from min point to max point)
    //               the choice between surf1 and surf2 for the envelope is
    //               the same for every point in the infinitesimal region
    //               above cv
    //               the surfaces are EQUAL over the curve cv
    Comparison_result
    operator()(const X_monotone_curve_2& cv,
               const Xy_monotone_surface_3& surf1,
               const Xy_monotone_surface_3& surf2) const {
      // a vertical surface cannot be defined in the infinitesimal region above
      // a curve
      CGAL_precondition(! surf1.is_vertical());
      CGAL_precondition(! surf2.is_vertical());

      CGAL_precondition(m_traits.compare_z_at_xy_3_object()
                        (cv, surf1, surf2) == EQUAL);
      CGAL_precondition(m_traits.compare_z_at_xy_3_object()
                        (cv.source(), surf1, surf2) == EQUAL);
      CGAL_precondition(m_traits.compare_z_at_xy_3_object()
                        (cv.target(), surf1, surf2) == EQUAL);

      if (m_traits.do_overlap(surf1, surf2)) return EQUAL;

      // now we must have 2 different non-vertical planes:
      // plane1: a1*x + b1*y + c1*z + d1 = 0  , c1 != 0
      // plane2: a2*x + b2*y + c2*z + d2 = 0  , c2 != 0

      const Plane_3& plane1 = surf1.plane();
      const Plane_3& plane2 = surf2.plane();

      FT a1 = plane1.a(), b1 = plane1.b(), c1 = plane1.c();
      FT a2 = plane2.a(), b2 = plane2.b(), c2 = plane2.c();

      // our line is a3*x + b3*y + c3 = 0
      // it is assumed that the planes intersect over this line
      const Line_2& line = cv.line();
      FT a3 = line.a(), b3 = line.b(), c3 = line.c();

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

      const Point_2& p1 = cv.left();
      const Point_2& p2 = cv.right();
      FT x1 = p1.x(), y1 = p1.y(), x2 = p2.x(), y2 = p2.y();

      Sign s2 = CGAL_NTS sign(-b3*x1+a3*y1-(-b3*x2+a3*y2));
      return s1 * s2;
    }
  };

  /*! Obtain a Compare_z_at_xy_above_3 functor object. */
  Compare_z_at_xy_above_3
  compare_z_at_xy_above_3_object() const
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
    using Traits_3 = Env_triangle_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Compare_z_at_xy_below_3(const Traits_3& traits) : m_traits(traits) {}

    friend class Env_triangle_traits_3<Kernel>;

  public:
    //
    Comparison_result
    operator()(const X_monotone_curve_2& cv,
               const Xy_monotone_surface_3& surf1,
               const Xy_monotone_surface_3& surf2) const {
      Comparison_result left_res =
        m_traits.compare_z_at_xy_above_3_object()(cv, surf1, surf2);
      return CGAL::opposite(left_res);

      /*if (left_res == LARGER)
        return SMALLER;
      else if (left_res == SMALLER)
        return LARGER;
      else
        return EQUAL;*/
    }
  };

  /*! Obtain a Compare_z_at_xy_below_3 functor object. */
  Compare_z_at_xy_below_3
  compare_z_at_xy_below_3_object() const
  { return Compare_z_at_xy_below_3(*this); }

  /***************************************************************************/

//  // checks if xy-monotone surface is vertical
//  class Is_vertical_3
//  {
//  public:
//
//    bool operator()(const Xy_monotone_surface_3& s) const
//    {
//      return false;
//    }
//  };
//
//  /*! Get a Is_vertical_3 functor object. */
//  Is_vertical_3 is_vertical_3_object() const
//  {
//    return Is_vertical_3();
//  }

  /***************************************************************************/

  // public method needed for testing

  // checks if point is in the xy-range of surf
  class Is_defined_over {
  protected:
    using Traits_3 = Env_triangle_traits_3<Kernel>;

    //! The traits (in case it has state).
    const Traits_3& m_traits;

    /*! Constructor
     * \param traits the traits
     */
    Is_defined_over(const Traits_3& traits) : m_traits(traits) {}

    friend class Env_triangle_traits_3<Kernel>;

  public:
    // checks if point is in the xy-range of surf
    bool operator()(const Point_2& point, const Xy_monotone_surface_3& surf)
      const {
      const Kernel& k = m_traits;

      // project the surface on the plane
      Triangle_2 boundary = m_traits.project(surf);

      // if surface is not vertical (i.e. boundary is not degenerate)
      // check if the projected point is inside the projected boundary
      if (! k.is_degenerate_2_object()(boundary))
        return (! k.has_on_unbounded_side_2_object()(boundary, point));

      // if surface is vertical, we check if the point is collinear
      // with the projected vertices, and on one of the projected segments
      // of the boundary
      Point_2 v1 = k.construct_vertex_2_object()(boundary, 0);
      Point_2 v2 = k.construct_vertex_2_object()(boundary, 1);
      Point_2 v3 = k.construct_vertex_2_object()(boundary, 2);

      if (! k.collinear_2_object()(v1, v2, point))
        return false;

      // enough to check 2 edges, because the 3rd is part of their union
      return (k.collinear_are_ordered_along_line_2_object()(v1, point, v2) ||
              k.collinear_are_ordered_along_line_2_object()(v2, point, v3));

    }
  };

  /*! Get a Is_defined_over functor object. */
  Is_defined_over is_defined_over_object() const
  { return Is_defined_over(*this); }

  //
  Segment_2 project (const Segment_3& seg) const {
    using Construct_vertex_3 = typename Kernel::Construct_vertex_3;

    const Kernel& k = *this;
    Construct_vertex_3 vertex_on = k.construct_vertex_3_object();
    const Point_3 q0 = vertex_on(seg, 0);
    const Point_3 q1 = vertex_on(seg, 1);
    const Point_2 p0(q0.x(), q0.y());
    const Point_2 p1(q1.x(), q1.y());
    return (k.construct_segment_2_object() (p0, p1));
  }

  //
  Point_2 project(const Point_3& obj) const
  { return Point_2(obj.x(), obj.y()); }

  //
  Triangle_2 project(const Xy_monotone_surface_3& triangle_3) const {
    const Point_3& end1 = triangle_3.vertex(0);
    const Point_3& end2 = triangle_3.vertex(1);
    const Point_3& end3 = triangle_3.vertex(2);
    Point_2 projected_end1(end1.x(), end1.y());
    Point_2 projected_end2(end2.x(), end2.y());
    Point_2 projected_end3(end3.x(), end3.y());
    return Triangle_2(projected_end1, projected_end2, projected_end3);
  }

  // triangles overlap if they lie on the same plane and intersect on it.
  // this test is only needed for non-vertical triangles
  bool do_overlap(const Xy_monotone_surface_3& s1,
                  const Xy_monotone_surface_3& s2) const {
    CGAL_precondition(s1.is_xy_monotone() && ! s1.is_vertical());
    CGAL_precondition(s2.is_xy_monotone() && ! s2.is_vertical());

    const Kernel& k = *this;
    auto do_x = k.do_intersect_3_object();
    if (! do_x(static_cast<Triangle_3>(s1), static_cast<Triangle_3>(s2)))
      return false;

    // check if they are coplanar
    Point_3 a1 = s1.vertex(0);
    Point_3 b1 = s1.vertex(1);
    Point_3 c1 = s1.vertex(2);
    Point_3 a2 = s2.vertex(0);
    Point_3 b2 = s2.vertex(1);
    Point_3 c2 = s2.vertex(2);
    bool b = k.coplanar_3_object()(a1, b1, c1, a2);
    if (! b) return false;

    b = k.coplanar_3_object()(a1, b1, c1, b2);
    if (! b) return false;

    b = k.coplanar_3_object()(a1, b1, c1, c2);
    return b;
  }

  // check whether two xy-monotone surfaces (3D-triangles or segments)
  // intersect
  bool do_intersect(const Xy_monotone_surface_3& s1,
                    const Xy_monotone_surface_3& s2) const {
    CGAL_precondition(s1.is_xy_monotone());
    CGAL_precondition(s2.is_xy_monotone());

    const Kernel& k = *this;
    auto do_x = k.do_intersect_3_object();
    if (! s1.is_segment() && ! s2.is_segment())
      return do_x(static_cast<Triangle_3>(s1), static_cast<Triangle_3>(s2));
    else if (! s1.is_segment())
      return do_x(static_cast<Triangle_3>(s1), static_cast<Segment_3>(s2));
    else if (! s2.is_segment())
      return do_x(static_cast<Segment_3>(s1), static_cast<Triangle_3>(s2));
    else return true; // if two segments, we don't use easy do-intersect test
  }

  // intersect two xy-monotone surfaces (3D-triangles or segments)
  // if the triangles overlap, the result is empty
  // the result can be a segment or a point
  std::optional<std::variant<Point_3, Segment_3>>
  intersection(const Xy_monotone_surface_3& s1, const Xy_monotone_surface_3& s2)
    const {
    CGAL_precondition(s1.is_xy_monotone());
    CGAL_precondition(s2.is_xy_monotone());
    Kernel k;

    // first, try to intersect the bounding boxes of the triangles,
    // efficiently return empty object when the triangles are faraway
    if (! CGAL::do_overlap(s1.bbox(), s2.bbox())) return std::nullopt;

    // if intersecting two segment - alculate the intersection
    // as in the case of dimension 2
    if (s1.is_segment() && s2.is_segment())
      return intersection_of_segments(s1, s2);

    // if both triangles lie on the same (non-vertical) plane, they overlap
    // we don't care about overlaps, because they are not passed to the
    // algorithm anyway, so we save the costly computation
    Plane_3 p1 = s1.plane();
    Plane_3 p2 = s2.plane();
    if ((p1 == p2) || (p1 == p2.opposite())) return std::nullopt;

    // calculate intersection between a triangle and the other triangle's
    // supporting plane
    // if there is no intersection - then the triangles have no intersection
    // between them.
    auto inter_obj = intersection(p1, s2);
    if (inter_obj == std::nullopt) return std::nullopt;

    // otherwise, if the intersection in a point, we should check if it lies
    // inside the first triangle
    if (const Point_3* inter_point = std::get_if<Point_3>(&(*inter_obj))) {
      std::optional<Point_3> res = intersection_on_plane_3(p1, s1, *inter_point);
      if (res != std::nullopt) return res.value();
    }
    else {
      // if the intersection is a segment, we check the intersection of the
      // other plane-triangle pair
      const Segment_3* inter_seg = std::get_if<Segment_3>(&(*inter_obj));
      CGAL_assertion(inter_seg != nullptr);

      auto inter_obj2 = intersection(p2, s1);

      // if there is no intersection - then the triangles have no intersection
      // between them.
      if (inter_obj2 == std::nullopt) return std::nullopt;

      if (const Point_3* inter_point = std::get_if<Point_3>(&(*inter_obj2))) {
        // if the intersection is a point, which lies on the segment,
        // than it is the result,
        // otherwise, empty result
         if (k.has_on_3_object()(*inter_seg, *inter_point)) return *inter_point;
         else return std::nullopt;
      }
      else {
        // both plane-triangle intersections are segments, which are collinear,
        // and lie on the line which is the intersection of the two supporting
        // planes
        const Segment_3* inter_seg2 = std::get_if<Segment_3>(&(*inter_obj));
        CGAL_assertion(inter_seg2 != nullptr);

        Point_3 min1 = k.construct_min_vertex_3_object()(*inter_seg);
        Point_3 max1 = k.construct_max_vertex_3_object()(*inter_seg);
        Point_3 min2 = k.construct_min_vertex_3_object()(*inter_seg2);
        Point_3 max2 = k.construct_max_vertex_3_object()(*inter_seg2);

         CGAL_assertion((k.collinear_3_object()(min1, min2, max1) &&
                         k.collinear_3_object()(min1, max2, max1)));

         // we need to find the overlapping part, if exists
         Point_3 min, max;
         if (k.less_xyz_3_object()(min1, min2)) min = min2;
         else min = min1;
         if (k.less_xyz_3_object()(max1, max2)) max = max1;
         else max = max2;

         Comparison_result comp_res = k.compare_xyz_3_object()(min, max);
         if (comp_res == EQUAL) return min;
         else if (comp_res == SMALLER) return Segment_3(min, max);
         // else - empty result
      }
    }
    return std::nullopt;
  }

  // calculate intersection between triangle & point on the same plane plane
  std::optional<Point_3>
  intersection_on_plane_3(const Plane_3& plane,
                          const Xy_monotone_surface_3& triangle,
                          const Point_3& point) const {
    const Kernel& k = *this;
    CGAL_precondition(triangle.is_xy_monotone() );
    CGAL_precondition(! k.is_degenerate_3_object()(plane) );
    CGAL_precondition((triangle.plane() == plane) ||
                      (triangle.plane() == plane.opposite()));
    CGAL_precondition(k.has_on_3_object()(plane, point) );
    CGAL_USE(plane);

    // if the point is inside the triangle, then the point is the intersection
    // otherwise there is no intersection
    bool has_on;
    if (triangle.is_segment())
      has_on = k.has_on_3_object()(static_cast<Segment_3>(triangle), point);
    else
      has_on = k.has_on_3_object()(static_cast<Triangle_3>(triangle), point);
    if (has_on) return point;
    else return std::nullopt;
  }

  // calculate intersection between 2 segments on the same vertical plane plane
  std::optional<std::variant<Point_3, Segment_3>>
  intersection_of_segments(const Xy_monotone_surface_3& s1,
                           const Xy_monotone_surface_3& s2) const {
    const Kernel& k = *this;
    CGAL_precondition(s1.is_xy_monotone() && s1.is_segment());
    CGAL_precondition(s2.is_xy_monotone() && s2.is_segment());

    // if the segments are not coplanar, they cannot intersect
    if (! k.coplanar_3_object()(s1.vertex(0), s1.vertex(1),
                                s2.vertex(0), s2.vertex(1)))
      return std::nullopt;

    const Plane_3& plane = s1.plane();
    if (s2.plane() != plane && s2.plane() != plane.opposite())
      // todo: this case is not needed in the algorithm,
      // so we don't implement it
      return std::nullopt;

    CGAL_precondition(! k.is_degenerate_3_object()(plane) );
    CGAL_precondition((s2.plane() == plane) ||
                      (s2.plane() == plane.opposite()));

    // for simplicity, we transform the segments to the xy-plane,
    // compute the intersection there, and transform it back to the 3d plane.
    Point_2 v1 = plane.to_2d(s1.vertex(0));
    Point_2 v2 = plane.to_2d(s1.vertex(1));
    Segment_2 seg1_t(v1, v2);

    Point_2 u1 = plane.to_2d(s2.vertex(0));
    Point_2 u2 = plane.to_2d(s2.vertex(1));
    Segment_2 seg2_t(u1, u2);

    auto inter_obj = k.intersect_2_object()(seg1_t, seg2_t);
    if (inter_obj == std::nullopt) return std::nullopt;

    if (const Point_2* inter_point = std::get_if<Point_2>(&(*inter_obj)))
      return plane.to_3d(*inter_point);

    const Segment_2* inter_segment = std::get_if<Segment_2>(&(*inter_obj));
    CGAL_assertion(inter_segment != nullptr);

    auto ctr_vertex = k.construct_vertex_2_object();
    return Segment_3(plane.to_3d(ctr_vertex(*inter_segment, 0)),
                     plane.to_3d(ctr_vertex(*inter_segment, 1)));
  }

  // calculate the intersection between a triangle/segment
  // and a (non degenerate) plane in 3d
  // the result object can be empty, a point, a segment or the original
  // triangle
  std::optional<std::variant<Point_3, Segment_3>>
  intersection(const Plane_3& pl, const Xy_monotone_surface_3& tri) const {
    const Kernel& k = *this;
    CGAL_precondition(tri.is_xy_monotone() );
    CGAL_precondition(! k.is_degenerate_3_object()(pl) );

    if (tri.is_segment())
      return k.intersect_3_object()(pl, static_cast<Segment_3>(tri));

    // first, check for all 3 vertices of tri on which side of pl they lie on
    int points_on_plane[3];    // contains the indices of vertices that lie
                               // on pl
    int points_on_positive[3]; // contains the indices of vertices that lie on
                               // the positive side of pl
    int points_on_negative[3]; // contains the indices of vertices that lie on
                               // the negative side of pl

    int n_points_on_plane = 0;
    int n_points_on_positive = 0;
    int n_points_on_negative = 0;

    Oriented_side side;
    for (int i = 0; i < 3; ++i) {
      side = pl.oriented_side(tri.vertex(i));
      if (side == ON_NEGATIVE_SIDE)
        points_on_negative[n_points_on_negative++] = i;
      else if (side == ON_POSITIVE_SIDE)
        points_on_positive[n_points_on_positive++] = i;
      else
        points_on_plane[n_points_on_plane++] = i;
    }

    CGAL_assertion(n_points_on_plane + n_points_on_positive +
                   n_points_on_negative == 3);

    // if all vertices of tri lie on the same size (positive/negative) of pl,
    // there is no intersection
    if ((n_points_on_positive == 3) || (n_points_on_negative == 3))
      return std::nullopt;

    // if all vertices of tri lie on pl then we return tri
    if (n_points_on_plane == 3) return tri;

    // if 2 vertices lie on pl, then return the segment between them
    if (n_points_on_plane == 2) {
      int point_idx1 = points_on_plane[0], point_idx2 = points_on_plane[1];
      return Segment_3(tri.vertex(point_idx1), tri.vertex(point_idx2));
    }

    // if only 1 lie on pl, should check the segment opposite to it on tri
    if (n_points_on_plane == 1) {
      int point_on_plane_idx = points_on_plane[0];

      // if the other 2 vertices are on the same side of pl,
      // then the answer is just this vertex
      if (n_points_on_negative == 2 || n_points_on_positive == 2)
        return tri.vertex(point_on_plane_idx);

      // now it is known that one vertex is on pl, and the segment of tri
      // opposite to it should intersect pl

      // the segment of tri opposite of tri[point_on_plane_idx]
      Segment_3 tri_segment(tri.vertex(point_on_plane_idx+1),
                            tri.vertex(point_on_plane_idx+2));

      auto inter_result = k.intersect_3_object()(pl, tri_segment);
      const Point_3* inter_point = std::get_if<Point_3>(&(*inter_result));
      CGAL_assertion( inter_point != nullptr );

      // create the resulting segment
      // (between tri[point_on_plane_idx] and inter_point)
      return Segment_3(tri.vertex(point_on_plane_idx), *inter_point);
    }

    CGAL_assertion(n_points_on_plane == 0);
    CGAL_assertion(n_points_on_positive + n_points_on_negative == 3);
    CGAL_assertion(n_points_on_positive != 0);
    CGAL_assertion(n_points_on_negative != 0);

    // now it known that there is an intersection between 2 segments of tri
    // and pl, it is also known which segments are those.
    Point_3 inter_points[2];
    int n_inter_points = 0;
    for (int pos_it = 0; pos_it < n_points_on_positive; ++pos_it)
      for (int neg_it = 0; neg_it < n_points_on_negative; ++neg_it) {
        Segment_3 seg(tri.vertex(points_on_positive[pos_it]),
                      tri.vertex(points_on_negative[neg_it]));
        auto inter_result = k.intersect_3_object()(pl, seg);
        const Point_3* inter_point = std::get_if<Point_3>(&(*inter_result));
        CGAL_assertion( inter_point != nullptr );
        inter_points[n_inter_points++] = *inter_point;
      }

    CGAL_assertion( n_inter_points == 2 );
    return Segment_3(inter_points[0], inter_points[1]);
  }

  // compare the value of s1 in p1 to the value of s2 in p2
  Comparison_result
  compare_z(const Point_2& p1, const Xy_monotone_surface_3& s1,
            const Point_2& p2, const Xy_monotone_surface_3& s2) {
    CGAL_precondition(is_defined_over_object()(p1, s1));
    CGAL_precondition(is_defined_over_object()(p2, s2));

    Point_3 v1 = envelope_point_of_surface(p1, s1);
    Point_3 v2 = envelope_point_of_surface(p2, s2);
    Kernel k;
    return k.compare_z_3_object()(v1, v2);
  }

  // find the envelope point of the surface over the given point
  // precondition: the surface is defined in point
  Point_3
  envelope_point_of_surface(const Point_2& p, const Xy_monotone_surface_3& s)
    const {
    CGAL_precondition(s.is_xy_monotone());
    CGAL_precondition(is_defined_over_object()(p, s));

    Point_3 point(p.x(), p.y(), 0);

    // Compute the intersetion between the vertical line and the given surfaces
    if (s.is_segment()) return envelope_point_of_segment(point, s);
    else {
      // s is a non-vertical triangle
      CGAL_assertion(! s.is_vertical());

      // Construct a vertical line passing through point
      Kernel k;
      Direction_3 dir (0, 0, 1);
      Line_3      vl = k.construct_line_3_object() (point, dir);

      const Plane_3& plane = s.plane();
      auto res = k.intersect_3_object()(plane, vl);
      CGAL_assertion(res != std::nullopt);
      const Point_3* ip = std::get_if<Point_3>(&(*res));
      CGAL_assertion(ip != nullptr);

      return *ip;
    }
  }

  // find the envelope point of the surface over the given point
  // precondition: the surface is defined in point and is a segment
  Point_3 envelope_point_of_segment(const Point_3& point,
                                    const Xy_monotone_surface_3& s) const {
    Kernel k;
    CGAL_precondition(s.is_segment());
    CGAL_precondition(is_defined_over_object()(project(point), s));

    // this is the vertical plane through the segment
    const Plane_3& plane = s.plane();

    // Construct a vertical line passing through point
    Direction_3 dir (0, 0, 1);
    Line_3 vl = k.construct_line_3_object() (point, dir);
    // we need 2 points on this line, to be transformed to 2d,
    // and preserve the direction of the envelope
    Point_3 vl_point1 = k.construct_point_on_3_object()(vl, 0);
    Point_3 vl_point2 = k.construct_point_on_3_object()(vl, 1);

    // the surface and the line are on the same plane(plane),
    // so we transform them to the xy-plane, compute the intersecting point
    // and transform it back to plane.
    const Point_3& v1 = s.vertex(0);
    const Point_3& v2 = s.vertex(1);

    Point_2 t1 = plane.to_2d(v1);
    Point_2 t2 = plane.to_2d(v2);

    Point_2 tvl_point1 = plane.to_2d(vl_point1);
    Point_2 tvl_point2 = plane.to_2d(vl_point2);
    Line_2 l(tvl_point1, tvl_point2);

    Segment_2 seg(t1, t2);
    auto inter_obj = k.intersect_2_object()(seg, l);
    const Point_2* inter_point = std::get_if<Point_2>(&(*inter_obj));
    CGAL_assertion(inter_point != nullptr);
    return plane.to_3d(*inter_point);
  }

  Point_2 construct_middle_point(const Point_2& p1, const Point_2& p2) const {
    Kernel k;
    return k.construct_midpoint_2_object()(p1, p2);
  }

  Point_2 construct_middle_point(const X_monotone_curve_2& cv) const {
    Kernel k;
    return k.construct_midpoint_2_object()(cv.source(), cv.target());
  }

  /***************************************************************************/
  // for vertical decomposition
  /***************************************************************************/

  class Construct_vertical_2 {
  public:
    X_monotone_curve_2 operator()(const Point_2& p1, const Point_2& p2) const
    { return X_monotone_curve_2(p1, p2); }
  };

  /*! Get a Construct_vertical_2 functor object. */
  Construct_vertical_2 construct_vertical_2_object() const
  { return Construct_vertical_2(); }

  Point_2 vertical_ray_shoot_2(const Point_2& pt, const X_monotone_curve_2& cv)
    const {
    CGAL_precondition(! cv.is_vertical());

    typename Kernel::Segment_2 seg = cv;
    const Kernel& k = *this;
    // If the curve contains pt, return it.
    if (k.has_on_2_object()(seg, pt)) return (pt);

    // Construct a vertical line passing through pt.
    typename Kernel::Direction_2 dir(0, 1);
    typename Kernel::Line_2 vl = k.construct_line_2_object()(pt, dir);

    // Compute the intersetion between the vertical line and the given curve.
    auto res = k.intersect_2_object()(seg, vl);
    const Point_2* ip = std::get_if<Point_2>(&(*res));
    CGAL_assertion(ip != nullptr);

    return *ip;
  }
};

/*! \class A representation of a triangle, as used by the
 * Env_triangle_traits_3 traits-class.
 */
template <typename Kernel_>
class Env_triangle_3 : public Env_triangle_traits_3<Kernel_>::_Triangle_cached_3
{
  using Kernel = Kernel_;
  using Triangle_3 = typename Kernel::Triangle_3;
  using Point_3 = typename Kernel::Point_3;
  using Plane_3 = typename Kernel::Plane_3;
  using Segment_3 = typename Kernel::Segment_3;

  using Base = typename Env_triangle_traits_3<Kernel>::_Triangle_cached_3;

public:
  /*! Default constructor.
   */
  Env_triangle_3() : Base() {}

  /*! Constructor from a "kernel" triangle.
   * \param seg The segment.
   */
  Env_triangle_3(const Triangle_3& tri) : Base(tri) {}

  /*! Construct a triangle from 3 end-points.
   * \param p1 The first point.
   * \param p2 The second point.
   * \param p3 The third point.
   */
  Env_triangle_3(const Point_3& p1, const Point_3& p2, const Point_3& p3) :
    Base(p1, p2, p3)
  {}

  /*! Construct a triangle from a plane and 3 end-points.
   * \param pl The supporting plane.
   * \param p1 The first point.
   * \param p2 The second point.
   * \param p3 The third point.
   * \pre All points must be on the supporting plane.
   */
  Env_triangle_3(const Plane_3& pl, const Point_3& p1,
                 const Point_3& p2, const Point_3& p3) :
    Base(pl, p1, p2, p3)
  {}

  /*! Construct a segment from 2 end-points.
   * \param p1 The first point.
   * \param p2 The second point.
   */
  Env_triangle_3(const Point_3& p1, const Point_3& p2) : Base(p1, p2) {}

  /*! Cast to a triangle.
   */
  operator Triangle_3() const
  { return (Triangle_3(this->vertex(0), this->vertex(1), this->vertex(2))); }

  /*! Cast to a segment (only when possible).
   */
  operator Segment_3() const {
    CGAL_precondition(this->is_segment());
    return (Segment_3(this->vertex(0), this->vertex(1)));
  }

  /*! Create a bounding box for the triangle.
   */
  Bbox_3 bbox() const {
    Triangle_3 tri(this->vertex(0), this->vertex(1), this->vertex(2));
    return (tri.bbox());
  }
};

template <typename Kernel_>
bool operator<(const Env_triangle_3<Kernel_>& a,
               const Env_triangle_3<Kernel_>& b) {
  if (a.vertex(0) < b.vertex(0)) return true;
  if (a.vertex(0) > b.vertex(0)) return false;
  if (a.vertex(1) < b.vertex(1)) return true;
  if (a.vertex(1) > b.vertex(1)) return false;
  if (a.vertex(2) < b.vertex(2)) return true;
  if (a.vertex(2) > b.vertex(2)) return false;

  return false;
}

template <typename Kernel_>
bool operator==(const Env_triangle_3<Kernel_>& a,
                const Env_triangle_3<Kernel_>& b) {
  return ((a.vertex(0) == b.vertex(0)) &&
          (a.vertex(1) == b.vertex(1)) &&
          (a.vertex(2) == b.vertex(2)));
}

/*! Exporter for the triangle class used by the traits-class.
 */
template <typename Kernel_, typename OutputStream>
OutputStream& operator<<(OutputStream& os, const Env_triangle_3<Kernel_>& tri) {
  os << static_cast<typename Kernel_::Triangle_3>(tri);
  if (tri.is_segment()) os << " (segment)";
  return os;
}

/*! Importer for the triangle class used by the traits-class.
 */
template <typename Kernel_, typename InputStream>
InputStream& operator>>(InputStream& is, Env_triangle_3<Kernel_>& tri) {
  typename Kernel_::Triangle_3 kernel_tri;
  is >> kernel_tri;
  tri = kernel_tri;
  return is;
}

} //namespace CGAL

#endif
