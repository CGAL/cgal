// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>
//                 Baruch Zukerman        <baruchzu@post.tau.ac.il>

/*! \file CGAL/Envelope_triangles_traits_3.h
 * \brief Model for CGAL's EnvelopeTraits_3 concept.
 * \endlink
 */

#ifndef CGAL_ENV_TRIANGLE_TRAITS_3_H
#define CGAL_ENV_TRIANGLE_TRAITS_3_H

#include <CGAL/license/Envelope_3.h>


#include <CGAL/Object.h>
#include <CGAL/enum.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Envelope_3/Envelope_base.h>

#include <vector>

namespace CGAL {

template <class Kernel_> class Env_triangle_3;

// this traits class supports both triagles and segments in 3d
template <class Kernel_>
class Env_triangle_traits_3 : public Arr_segment_traits_2<Kernel_>
{
public:
  typedef Arr_segment_traits_2<Kernel_>             Traits_2;
  typedef typename Traits_2::Point_2                Point_2;
  typedef typename Traits_2::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Traits_2::Multiplicity           Multiplicity;

  typedef Kernel_                                   Kernel;
  typedef Env_triangle_traits_3<Kernel>             Self;

  typedef typename Kernel::Point_3                  Point_3;

  /*!
   * \class Representation of a 3d triangle with cached data.
   */
  class _Triangle_cached_3 
  {
  public:

    typedef typename Kernel::Plane_3               Plane_3;
    typedef typename Kernel::Triangle_3            Triangle_3;
    typedef typename Kernel::Point_3               Point_3;
    typedef typename Kernel::Segment_3             Segment_3;

  protected:

    Plane_3 pl;          // The plane that supports the triangle.
    Point_3 vertices[3]; // The vertices of the triangle.
    bool    is_vert;     // Is this a vertical triangle (or a segment).
    bool    is_seg;  // Is this a segment.
  public:

    /*!
     * Default constructor.
     */
    _Triangle_cached_3() :
      is_vert(false),
      is_seg(false)
    {}

    /*!
     * Constructor from a non-degenerate triangle.
     * \param tri The triangle.
     * \pre The triangle is not degenerate.
     */
    _Triangle_cached_3(const Triangle_3 & tri)
    {
      Kernel   kernel;
      CGAL_assertion(!kernel.is_degenerate_3_object()(tri));

      typename Kernel::Construct_vertex_3
        construct_vertex = kernel.construct_vertex_3_object();

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

    /*!
     * Construct a triangle from three non-collinear end-points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \param p3 The third point.
     * \pre The 3 endpoints are not the collinear.
     */
    _Triangle_cached_3(const Point_3 &p1, const Point_3 &p2,
                       const Point_3 &p3)
    {
      Kernel   kernel;
      CGAL_assertion(!kernel.collinear_3_object()(p1, p2, p3));
      
      vertices[0] = p1;
      vertices[1] = p2;
      vertices[2] = p3;

      pl = kernel.construct_plane_3_object()(vertices[0],
                                             vertices[1],
                                             vertices[2]);
      Self self;
      is_vert = kernel.collinear_2_object()(self.project(vertices[0]),
                                            self.project(vertices[1]),
                                            self.project(vertices[2]));
      is_seg = false;
    }

    /*!
     * Construct a triangle from 3 end-points on a supporting plane.
     * \param supp_plane The supporting plane.
     * \param p1 The first point.
     * \param p2 The second point.
     * \param p3 The third point.
     * \pre The 3 endpoints are not the collinear and all lie on the given
     *      plane.
     */
    _Triangle_cached_3(const Plane_3& supp_plane,
                       const Point_3 &p1,
                       const Point_3 &p2,
                       const Point_3 &p3) :
      pl(supp_plane)
    {
      Kernel   kernel;

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

    /*!
     * Constructor from a segment.
     * \param seg The segment.
     * \pre The segment is not degenerate.
     */
    _Triangle_cached_3(const Segment_3 & seg)
    {
      Kernel   kernel;
      CGAL_assertion(!kernel.is_degenerate_3_object()(seg));

      typename Kernel::Construct_vertex_3
        construct_vertex = kernel.construct_vertex_3_object();

      vertices[0] = construct_vertex(seg, 0);
      vertices[1] = construct_vertex(seg, 1);
      vertices[2] = vertices[1];

      is_vert = true;
      is_seg = true;

      // construct a vertical plane through the segment
      Point_3 tmp(vertices[0].x(), vertices[0].y(), vertices[0].z()-1);
      pl = kernel.construct_plane_3_object()(vertices[0],
                                             vertices[1], tmp);

    }

    /*!
     * Constructor from two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \param seg The segment.
     * \pre The segment between the points is not degenerate.
     */
    _Triangle_cached_3(const Point_3 &p1, const Point_3 &p2)
    {
      Kernel   kernel;
      CGAL_assertion(!kernel.equal_3_object()(p1, p2));

      vertices[0] = p1;
      vertices[1] = p2;
      vertices[2] = p2;

      is_vert = true;
      is_seg = true;

      // construct a vertical plane through the segment
      Point_3 tmp(vertices[0].x(), vertices[0].y(), vertices[0].z()-1);
      pl = kernel.construct_plane_3_object()(vertices[0],
                                             vertices[1], tmp);

    }

    /*!
     * Assignment operator.
     * \param tri the source triangle to copy from
     */
    const _Triangle_cached_3& operator=(const Triangle_3 &tri)
    {
      Kernel   kernel;
      CGAL_assertion(!kernel.is_degenerate_3_object()(tri));

      typename Kernel_::Construct_vertex_3
        construct_vertex = kernel.construct_vertex_3_object();

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

    /*!
     * Get the ith endpoint.
     */
    const Point_3& vertex(unsigned int i) const
    {
      return vertices[i%3];
    }

    /*!
     * Get the supporting plane.
     */
    const Plane_3& plane() const
    {
      return (pl);
    }

    /*!
     * Check if the triangel is vertical.
     */
    bool is_vertical() const
    {
      return (is_vert);
    }

    /*!
     * Check if the surface is a segment.
     */
    bool is_segment() const
    {
      return (is_seg);
    }

    /*!
     * Check if the surface is xy-monotone (false, if it is a vertical
     * triangle)
     */
    bool is_xy_monotone() const
    {
      return (!is_vertical() || is_segment());
    }
  };

public:

  // types for EnvelopeTraits_3 concept
  //! type of xy-monotone surfaces
  typedef Env_triangle_3<Kernel>                    Xy_monotone_surface_3;
  //! type of surfaces
  typedef Xy_monotone_surface_3                     Surface_3;

  // we have a collision between the Kernel's Intersect_2 and the one
  // from the segment traits
  typedef typename Traits_2::Intersect_2            Intersect_2;

protected:
  typedef typename Kernel::FT                       FT;
  typedef typename Kernel::Triangle_2               Triangle_2;
  typedef typename Kernel::Segment_2                Segment_2;

  typedef typename Kernel::Segment_3                Segment_3;
  typedef typename Kernel::Triangle_3               Triangle_3;
  typedef typename Kernel::Plane_3                  Plane_3;

  typedef typename Kernel::Assign_2                 Assign_2;
  typedef typename Kernel::Construct_vertex_2       Construct_vertex_2;

  typedef typename Kernel::Assign_3                 Assign_3;
  typedef typename Kernel::Intersect_3              Intersect_3;
  typedef typename Kernel::Construct_vertex_3       Construct_vertex_3;


  typedef typename Kernel::Line_2                   Line_2;
  typedef typename Kernel::Direction_2              Direction_2;

  typedef typename Kernel::Line_3                   Line_3;
  typedef typename Kernel::Direction_3              Direction_3;

  typedef std::pair<X_monotone_curve_2,
                    Multiplicity>                   Intersection_curve;
public:

  /***************************************************************************/
  // EnvelopeTraits_3 functors
  /***************************************************************************/

  /*!\brief
   * Subdivide the given surface into envelope relevant xy-monotone 
   * parts, and insert them into the output iterator.
   * 
   * The iterator value-type is Xy_monotone_surface_3
   */
  class Make_xy_monotone_3
  {
  protected:
    const Self * parent;

  public:
    Make_xy_monotone_3(const Self * p) : parent(p)
    {}
    // create xy-monotone surfaces from a general surface
    // return a past-the-end iterator
    template <class OutputIterator>
    OutputIterator operator()(const Surface_3& s,
                              bool is_lower,
                              OutputIterator o) const
    {
      m_is_lower = is_lower;

      // a non-vertical triangle is already xy-monotone
      if (!s.is_vertical())
        *o++ = s;
      else
      {        
        // split a vertical triangle into one or two segments
        const Point_3 &a1 = s.vertex(0),
                       a2 = s.vertex(1),
                       a3 = s.vertex(2);
        Point_2 b1 = parent->project(a1),
                b2 = parent->project(a2),
                b3 = parent->project(a3);
        Kernel k;
        if (k.collinear_are_ordered_along_line_2_object()(b1, b2, b3))
        {
          if (k.equal_2_object()(b1, b2))
            // only one segment in the output - the vertical does not count
            *o++ = Xy_monotone_surface_3(find_envelope_point(a1, a2), a3);
          else if (k.equal_2_object()(b2, b3))
            *o++ = Xy_monotone_surface_3(a1, find_envelope_point(a2, a3));
          else
            // check whether two or one segments appear on the envelope
            return find_envelope_segments(a1, a2, a3, s.plane(), o);
        }
        else if (k.collinear_are_ordered_along_line_2_object()(b1, b3, b2))
        {
          if (k.equal_2_object()(b1, b3))
            // only one segment in the output
            *o++ = Xy_monotone_surface_3(find_envelope_point(a1, a3), a2);
          else
            // check whether two or one segments appear on the envelope
            return find_envelope_segments(a1, a3, a2, s.plane(), o);
        }
        else
        {
          // check whether two or one segments appear on the envelope
          return find_envelope_segments(a2, a1, a3, s.plane(), o);
        }

      }
      return o;
    }

  protected:
    // find the envelope point among the two points with same xy coordinates
    const Point_3& find_envelope_point (const Point_3& p1,
                                        const Point_3& p2) const
    {
      CGAL_precondition(p1.x() == p2.x() && p1.y() == p2.y());
      Kernel k;
      Comparison_result cr = k.compare_z_3_object()(p1, p2);
      CGAL_assertion(cr != EQUAL);
      if ((m_is_lower && cr == SMALLER) ||
          (!m_is_lower && cr == LARGER))
        return p1;
      else
        return p2;      
    }

    // get the three triangle coordinates (ordered along 2d-line) and find
    // which segment(s) is(are) the envelope of this triangle
    // "plane" is the vertical plane on which the triangle lies
    template <class OutputIterator>
    OutputIterator find_envelope_segments(const Point_3& p1,
                                          const Point_3& p2,
                                          const Point_3& p3,
                                          const Plane_3& plane,
                                          OutputIterator o) const
    {
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
      FT w1 = a*p1.y() - b*p1.x(),
         w2 = a*p2.y() - b*p2.x(),
         w3 = a*p3.y() - b*p3.x();
      
      Sign s1 = CGAL::sign((w3 - w1)*(p2.z() - p1.z()) - 
                           (p3.z() - p1.z())*(w2 - w1));

      // the points should not be collinear
      CGAL_assertion(s1 != 0);

      // should also take care for the original and trasformed direction of
      // the segment
      Sign s2 = CGAL_NTS sign(w3 - w1);
      Sign s = CGAL_NTS sign(int(s1 * s2));
                  
      bool use_one_segment = true;
      if ((m_is_lower  && s == NEGATIVE) ||
          (!m_is_lower && s == POSITIVE))
        use_one_segment = false;

      if (use_one_segment)
      {
        *o++ = Xy_monotone_surface_3(p1, p3);
      }
      else
      {
        *o++ = Xy_monotone_surface_3(p1, p2);
        *o++ = Xy_monotone_surface_3(p2, p3);
      }
      return o;
    }

    mutable bool m_is_lower;
  };

  /*! Get a Make_xy_monotone_3 functor object. */
  Make_xy_monotone_3
  make_xy_monotone_3_object() const
  {
    return Make_xy_monotone_3(this);
  }

  /*!\brief
   * Insert all 2D curves, which form the boundary of the vertical
   * projection of the surface onto the xy-plane, into the output iterator.
   * The iterator value-type is X_monotone_curve_2.
   */
  class Construct_projected_boundary_2
  {
  protected:
    const Self *parent;
  public:

    Construct_projected_boundary_2(const Self* p)
      : parent(p)
    {}

    // insert into the OutputIterator all the (2d) curves of the boundary of
    // the vertical projection of the surface on the xy-plane
    // the OutputIterator value type is X_monotone_curve_2
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s,
                              OutputIterator o) const
    {
      // the input xy-monotone surface should be either non-vertical or
      // a segment
      CGAL_assertion(s.is_xy_monotone());
      
      if (!s.is_vertical())
      {
        // the projection is a triangle
        const Point_3 &a1 = s.vertex(0),
                       a2 = s.vertex(1),
                       a3 = s.vertex(2);
        Point_2 b1 = parent->project(a1),
                b2 = parent->project(a2),
                b3 = parent->project(a3);

        Kernel k;

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

        CGAL_assertion(s1 != ON_ORIENTED_BOUNDARY && 
                       s2 != ON_ORIENTED_BOUNDARY &&
                       s3 != ON_ORIENTED_BOUNDARY);
        
        *o++ = make_object(std::make_pair(A, s1));
        *o++ = make_object(std::make_pair(B, s2));
        *o++ = make_object(std::make_pair(C, s3));
      }
      else
      {
        // s is a segment, and so is its projection
        // s shouldn't be a z-vertical segment
        const Point_3 &a1 = s.vertex(0),
                       a2 = s.vertex(1);
        
        Point_2 b1 = parent->project(a1),
                b2 = parent->project(a2);
        CGAL_assertion(b1 != b2);
                
        *o++ = make_object(std::make_pair(X_monotone_curve_2(b1, b2), 
                                          ON_ORIENTED_BOUNDARY));
      }
      return o;
    }  
  };  
  
  /*! Get a Construct_projected_boundary_curves_2 functor object. */
  Construct_projected_boundary_2
  construct_projected_boundary_2_object() const
  {
    return Construct_projected_boundary_2(this);
  }

  /*!\brief
   * Insert all the 2D projections (onto the xy-plane) of the 
   * intersection objects between s1 and s2 into the output iterator.
   *
   * The iterator value-type is Object. An Object may be:
   * 1. A pair<X_monotone_curve_2,Intersection_type>, where the intersection 
   * type is an enumeration that can take the values
   * {Transversal, Tangency, Unknown}.
   * 2. A Point_2 instance (in degenerate cases).
   */
  class Construct_projected_intersections_2
  {
  protected:
    const Self *parent;
  public:

    Construct_projected_intersections_2(const Self* p)
      : parent(p)
    {}
    
    // insert into OutputIterator all the (2d) projections on the xy plane of
    // the intersection objects between the 2 surfaces
    // the data type of OutputIterator is Object
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s1,
                              const Xy_monotone_surface_3& s2,
                              OutputIterator o) const
    {
      CGAL_assertion(s1.is_xy_monotone() && s2.is_xy_monotone());
      
      Kernel k;
      if (!parent->do_intersect(s1, s2))
      {
        return o;
      }
        
      Object inter_obj = parent->intersection(s1,s2);
      if (inter_obj.is_empty())
      {
        return o;
      }

      Point_3 point;
      Segment_3 curve;
      if (k.assign_3_object()(point, inter_obj))
        *o++ = make_object(parent->project(point));
      else
      {
        CGAL_assertion_code(bool b = )
        k.assign_3_object()(curve, inter_obj);
        CGAL_assertion(b);

        Segment_2  proj_seg = parent->project(curve);
        if (! k.is_degenerate_2_object() (proj_seg))
        {
          Intersection_curve inter_cv (proj_seg, 1);
          *o++ = make_object(inter_cv);
        }
        else
        {
          const Point_2&  p = k.construct_point_on_2_object() (proj_seg, 0);
          *o++ = make_object(p);
        }
      }

      return o;
    }  
  };  

  /*! Get a Construct_projected_intersections_2 functor object. */
  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const
  {
    return Construct_projected_intersections_2(this);
  }

  /*!\brief
   * Check if the surface s1 is closer/equally distanced/farther 
   * from the envelope with respect to s2 at the xy-coordinates of p/c.
   */
  class Compare_z_at_xy_3
  {
  protected:
    const Self *parent;

  public:

    Compare_z_at_xy_3(const Self* p)
      : parent(p)
    {}

    // check which of the surfaces is closer to the envelope at the xy 
    // coordinates of point
    // (i.e. lower if computing the lower envelope, or upper if computing 
    // the upper envelope)
    // precondition: the surfaces are defined in point
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& surf1,
                                 const Xy_monotone_surface_3& surf2) const
    {      
      // we compute the points on the planes, and then compare their z 
      // coordinates
      const Plane_3& plane1 = surf1.plane();
      const Plane_3& plane2 = surf2.plane();

      // if the 2 triangles have the same supporting plane, and they are not 
      // vertical, then they have the same z coordinate over this point
      if ((plane1 == plane2 || plane1 == plane2.opposite()) &&
          !surf1.is_vertical())
      {
        return EQUAL;
      }

      Kernel k;

      // Compute the intersetion between the vertical line and the given 
      // surfaces
      Point_3 ip1 = parent->envelope_point_of_surface(p, surf1);
      Point_3 ip2 = parent->envelope_point_of_surface(p, surf2);
      
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
                                 const Xy_monotone_surface_3& surf2) const
    {      
      // first try the endpoints, if cannot be sure, use the mid point
      Comparison_result res;
      res = parent->compare_z_at_xy_3_object()(cv.left(), surf1, surf2);

      if (res == EQUAL)
      {
        res = parent->compare_z_at_xy_3_object()(cv.right(), surf1, surf2);
        if (res == EQUAL)
        {
          Point_2 mid = parent->construct_middle_point(cv);
          res = parent->compare_z_at_xy_3_object()(mid, surf1, surf2);
        }
      }
      
      return res;
    }
  
  };
   
  /*! Get a Compare_z_at_xy_3 functor object. */
  Compare_z_at_xy_3 
  compare_z_at_xy_3_object() const
  {
    return Compare_z_at_xy_3(this);
  }

  /*!\brief 
   * Check if the surface s1 is closer/equally distanced/farther 
   * from the envelope with
   * respect to s2 immediately above the curve c. 
   */
  class Compare_z_at_xy_above_3
  {
  protected:
    const Self *parent;

  public:

    Compare_z_at_xy_above_3(const Self* p)
      : parent(p)
    {}
    
    // check which of the surfaces is closer to the envelope on the points 
    // above the curve cv
    // (i.e. lower if computing the lower envelope, or upper if computing the
    // upper envelope)
    // precondition: the surfaces are defined above cv (to the left of cv, 
    //               if cv is directed from min point to max point)
    //               the choise between surf1 and surf2 for the envelope is 
    //               the same for every point in the infinitesimal region 
    //               above cv 
    //               the surfaces are EQUAL over the curve cv
    Comparison_result
    operator()(const X_monotone_curve_2& cv,
               const Xy_monotone_surface_3& surf1,
               const Xy_monotone_surface_3& surf2) const
    {
      // a vertical surface cannot be defined in the infinitesimal region above
      // a curve
      CGAL_precondition(!surf1.is_vertical());
      CGAL_precondition(!surf2.is_vertical());

      CGAL_precondition(parent->compare_z_at_xy_3_object()
                              (cv, surf1, surf2) == EQUAL);
      CGAL_precondition(parent->compare_z_at_xy_3_object()
                              (cv.source(), surf1, surf2) == EQUAL);
      CGAL_precondition(parent->compare_z_at_xy_3_object()
                              (cv.target(), surf1, surf2) == EQUAL);

      
      if (parent->do_overlap(surf1, surf2))
      {
        return EQUAL;
      }

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


  /*! Get a Compare_z_at_xy_above_3 functor object. */
  Compare_z_at_xy_above_3
  compare_z_at_xy_above_3_object() const
  {
    return Compare_z_at_xy_above_3(this);
  }

  /*!\brief 
   * Check if the surface s1 is closer/equally distanced/farther 
   * from the envelope with
   * respect to s2 immediately below the curve c. 
   */
  class Compare_z_at_xy_below_3
  {
  protected:
    const Self *parent;

  public:

    Compare_z_at_xy_below_3(const Self* p)
      : parent(p)
    {}
    
    Comparison_result
    operator()(const X_monotone_curve_2& cv,
               const Xy_monotone_surface_3& surf1,
               const Xy_monotone_surface_3& surf2) const
    {
      Comparison_result left_res = 
        parent->compare_z_at_xy_above_3_object()(cv, surf1, surf2);
      return CGAL::opposite(left_res);

      /*if (left_res == LARGER)
        return SMALLER;
      else if (left_res == SMALLER)
        return LARGER;
      else
        return EQUAL;*/
    }  
  };

  /*! Get a Compare_z_at_xy_below_3 functor object. */
  Compare_z_at_xy_below_3
  compare_z_at_xy_below_3_object() const
  {
    return Compare_z_at_xy_below_3(this);
  }

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
  class Is_defined_over
  {
  public:
    // checks if point is in the xy-range of surf
    bool operator()(const Point_2& point, 
		    const Xy_monotone_surface_3& surf) const

    {
      Kernel k;
      Self parent;

      // project the surface on the plane
      Triangle_2 boundary = parent.project(surf);

      // if surface is not vertical (i.e. boundary is not degenerate)
      // check if the projected point is inside the projected boundary
      if (!k.is_degenerate_2_object()(boundary))
        return (!k.has_on_unbounded_side_2_object()(boundary, point));

      // if surface is vertical, we check if the point is collinear
      // with the projected vertices, and on one of the projected segments
      // of the boundary
      Point_2 v1 = k.construct_vertex_2_object()(boundary, 0);
      Point_2 v2 = k.construct_vertex_2_object()(boundary, 1);
      Point_2 v3 = k.construct_vertex_2_object()(boundary, 2);

      if (!k.collinear_2_object()(v1, v2, point))
        return false;

      // enough to check 2 edges, because the 3rd is part of their union
      return (k.collinear_are_ordered_along_line_2_object()(v1, point, v2) ||
              k.collinear_are_ordered_along_line_2_object()(v2, point, v3));

    }
  };

  /*! Get a Is_defined_over functor object. */
  Is_defined_over is_defined_over_object() const
  {
    return Is_defined_over();
  }

  Segment_2 project (const Segment_3& seg) const
  {
    typedef typename Kernel::Construct_vertex_3 Construct_vertex_3;
    
    Kernel              k;
    Construct_vertex_3  vertex_on = k.construct_vertex_3_object();

    const Point_3      q0 = (vertex_on (seg, 0));
    const Point_3      q1 = (vertex_on (seg, 1));
    const Point_2      p0 (q0.x(), q0.y());
    const Point_2      p1 (q1.x(), q1.y());
    
    return (k.construct_segment_2_object() (p0, p1));
  }

  Point_2 project(const Point_3& obj) const
  {
    return Point_2(obj.x(), obj.y());
  }
  
  Triangle_2 project(const Xy_monotone_surface_3& triangle_3) const
  {
    const Point_3&  end1 = triangle_3.vertex(0),
                    end2 = triangle_3.vertex(1),
                    end3 = triangle_3.vertex(2);
    Point_2 projected_end1(end1.x(), end1.y()),
            projected_end2(end2.x(), end2.y()),

            projected_end3(end3.x(), end3.y());
    return Triangle_2(projected_end1, projected_end2, projected_end3);
  }

  // triangles overlap if they lie on the same plane and intersect on it.
  // this test is only needed for non-vertical triangles
  bool do_overlap(const Xy_monotone_surface_3& s1, 
		              const Xy_monotone_surface_3& s2) const
  {
    CGAL_precondition(s1.is_xy_monotone() && !s1.is_vertical());
    CGAL_precondition(s2.is_xy_monotone() && !s2.is_vertical());

    Kernel k;
    if (!k.do_intersect_3_object()(static_cast<Triangle_3>(s1),
                                   static_cast<Triangle_3>(s2)))
      return false;

    // check if they are coplanar
    Point_3 a1 = s1.vertex(0),
            b1 = s1.vertex(1),
            c1 = s1.vertex(2);
    Point_3 a2 = s2.vertex(0),
            b2 = s2.vertex(1),
            c2 = s2.vertex(2);
    bool b = k.coplanar_3_object()(a1, b1, c1, a2);
    if (!b) return false;

    b = k.coplanar_3_object()(a1, b1, c1, b2);
    if (!b) return false;

    b = k.coplanar_3_object()(a1, b1, c1, c2);
    return b;    
  }

  // check whethe two xy-monotone surfaces (3D-triangles or segments)
  // intersect
  bool do_intersect(const Xy_monotone_surface_3& s1,
                    const Xy_monotone_surface_3& s2) const
  {
    CGAL_precondition(s1.is_xy_monotone());
    CGAL_precondition(s2.is_xy_monotone());
    Kernel k;
    if (!s1.is_segment() && !s2.is_segment())
      return k.do_intersect_3_object()(static_cast<Triangle_3>(s1),
                                       static_cast<Triangle_3>(s2));
    else if (!s1.is_segment())
      return k.do_intersect_3_object()(static_cast<Triangle_3>(s1),
                                       static_cast<Segment_3>(s2));
    else if (!s2.is_segment())
      return k.do_intersect_3_object()(static_cast<Segment_3>(s1),
                                       static_cast<Triangle_3>(s2));
    else
      // in case of two segments, we don't use easy do-intersect test
      return true;
  }
  
  // intersect two xy-monotone surfaces (3D-triangles or segments)
  // if the triangles overlap, the result is empty
  // the result can be a segment or a point
  Object intersection(const Xy_monotone_surface_3& s1, 
                      const Xy_monotone_surface_3& s2) const
  {
    CGAL_precondition(s1.is_xy_monotone());
    CGAL_precondition(s2.is_xy_monotone());
    Kernel k;

    // first, try to intersect the bounding boxes of the triangles,
    // efficiently return empty object when the triangles are faraway
    if (!CGAL::do_overlap(s1.bbox(), s2.bbox()))
      return Object();

    // if intersecting two segment - alculate the intersection
    // as in the case of dimention 2
    if (s1.is_segment() && s2.is_segment())
    {
      Object res = intersection_of_segments(s1, s2);
      return res;
    }
  
    // if both triangles lie on the same (non-vertical) plane, they overlap
    // we don't care about overlaps, because they are not passed to the
    // algorithm anyway, so we save the costly computation
    Plane_3 p1 = s1.plane();
    Plane_3 p2 = s2.plane();
    if  (p1 == p2 || p1 == p2.opposite())
        return Object();

    // calculate intersection between a triangle and the other triangle's 
    // supporting plane
    // if there is no intersection - then the triangles have no intersection 
    // between them.
    Object inter_obj = intersection(p1, s2);
      
    if (inter_obj.is_empty())
      return Object();

    // otherwise, if the intersection in a point, we should check if it lies
    // inside the first triangle
    Assign_3 assign_obj = k.assign_3_object();
    Point_3 inter_point;
    if (assign_obj(inter_point, inter_obj))
    {
      Object res = intersection_on_plane_3(p1, s1, inter_point);
      return res;
    }
    else
    {
      // if the intersection is a segment, we check the intersection of the
      // other plane-triangle pair
      Segment_3 inter_seg;
      CGAL_assertion(assign_obj(inter_seg, inter_obj));
      assign_obj(inter_seg, inter_obj);

      inter_obj = intersection(p2, s1);

      // if there is no intersection - then the triangles have no intersection 
      // between them.
      if (inter_obj.is_empty())
      	return Object();
      
      if (assign_obj(inter_point, inter_obj))
      {
      	// if the intersection is a point, which lies on the segment,
      	// than it is the result,
      	// otherwise, empty result
      	 if (k.has_on_3_object()(inter_seg, inter_point))
      	   return make_object(inter_point);
      	 else
      	   return Object();
      }
      else
      {
      	// both plane-triangle intersections are segments, which are collinear,
      	// and lie on the line which is the intersection of the two supporting
      	// planes
        Segment_3 inter_seg2;
      	CGAL_assertion(assign_obj(inter_seg2, inter_obj));
      	assign_obj(inter_seg2, inter_obj);
	
      	Point_3 min1 = k.construct_min_vertex_3_object()(inter_seg),
      	        max1 = k.construct_max_vertex_3_object()(inter_seg);
      	Point_3 min2 = k.construct_min_vertex_3_object()(inter_seg2),
      	        max2 = k.construct_max_vertex_3_object()(inter_seg2); 

       	CGAL_assertion((k.collinear_3_object()(min1, min2, max1) &&
                         k.collinear_3_object()(min1, max2, max1)));

       	// we need to find the overlapping part, if exists
       	Point_3 min, max;
       	if (k.less_xyz_3_object()(min1, min2))
       	  min = min2;
       	else
       	  min = min1;
       	if (k.less_xyz_3_object()(max1, max2))
       	  max = max1;
       	else
       	  max = max2;

       	Object res;
       	Comparison_result comp_res = k.compare_xyz_3_object()(min, max);
       	if (comp_res == EQUAL)
       	  res = make_object(min);
       	else if (comp_res == SMALLER)
       	  res = make_object(Segment_3(min, max));
       	// else - empty result

       	return res;
      }
    }
  }

  // calculate intersection between triangle & point on the same plane plane
  Object intersection_on_plane_3(const Plane_3& plane,
                                 const Xy_monotone_surface_3& triangle,
                                 const Point_3& point) const
  {
    Kernel k;
    CGAL_precondition( triangle.is_xy_monotone() );
    CGAL_precondition( !k.is_degenerate_3_object()(plane) );
    CGAL_precondition( triangle.plane() == plane ||
                       triangle.plane() == plane.opposite());
    CGAL_precondition( k.has_on_3_object()(plane, point) );
    CGAL_USE(plane);

    // if the point is inside the triangle, then the point is the intersection
    // otherwise there is no intersection
    bool has_on;
    if (triangle.is_segment())
      has_on = k.has_on_3_object()(static_cast<Segment_3>(triangle), point);
    else
      has_on = k.has_on_3_object()(static_cast<Triangle_3>(triangle), point);
    if (has_on)
      return make_object(point);
    else
      return Object();
  }

  // calculate intersection between 2 segments on the same vertical plane plane
  Object intersection_of_segments(const Xy_monotone_surface_3& s1,
                                  const Xy_monotone_surface_3& s2) const
  {
    Kernel k;
    CGAL_precondition( s1.is_xy_monotone() && s1.is_segment());
    CGAL_precondition( s2.is_xy_monotone() && s2.is_segment());

    // if the segments are not coplanar, they cannot intersect
    if (!k.coplanar_3_object()(s1.vertex(0), s1.vertex(1),
                               s2.vertex(0), s2.vertex(1)))
      return Object();

    const Plane_3& plane = s1.plane();
    if (s2.plane() != plane &&
        s2.plane() != plane.opposite())
      // todo: this case is not needed in the algorithm,
      // so we don't implement it
      return Object();

    CGAL_precondition( !k.is_degenerate_3_object()(plane) );
    CGAL_precondition( s2.plane() == plane ||
                       s2.plane() == plane.opposite());

    // for simplicity, we transform the segments to the xy-plane,
    // compute the intersection there, and transform it back to the 3d plane.
    Point_2 v1 = plane.to_2d(s1.vertex(0)),
            v2 = plane.to_2d(s1.vertex(1));
    Segment_2 seg1_t(v1, v2);

  	Point_2 u1 = plane.to_2d(s2.vertex(0)),
            u2 = plane.to_2d(s2.vertex(1));
  	Segment_2 seg2_t(u1, u2);

  	Object inter_obj = k.intersect_2_object()(seg1_t, seg2_t);
  	Assign_2 assign_2 = k.assign_2_object();
  	if (inter_obj.is_empty())
  		return inter_obj;

  	Point_2 inter_point;
    Segment_2 inter_segment;

    if (assign_2(inter_point, inter_obj))
  	  return make_object(plane.to_3d(inter_point));
    else
    {
      CGAL_assertion_code(bool b = )
      assign_2(inter_segment, inter_obj);
      CGAL_assertion(b);
      
      return make_object 
        (Segment_3
         (plane.to_3d(k.construct_vertex_2_object()(inter_segment, 0)),
          plane.to_3d(k.construct_vertex_2_object()(inter_segment, 1))));
    }

  }

  // calculate the intersection between a triangle/segment
  // and a (non degenerate) plane in 3d
  // the result object can be empty, a point, a segment or the original
  // triangle
  Object intersection(const Plane_3& pl, 
		                  const Xy_monotone_surface_3& tri) const
  {
    Kernel k;
    CGAL_precondition( tri.is_xy_monotone() );
    CGAL_precondition( !k.is_degenerate_3_object()(pl) );

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
    for (int i=0; i<3; ++i)
    {
      side = pl.oriented_side(tri.vertex(i));
      if (side == ON_NEGATIVE_SIDE)
        points_on_negative[n_points_on_negative++] = i;
      else if (side == ON_POSITIVE_SIDE)
        points_on_positive[n_points_on_positive++] = i;
      else
        points_on_plane[n_points_on_plane++] = i;
    }

    CGAL_assertion(n_points_on_plane + 
            n_points_on_positive + n_points_on_negative == 3);

    // if all vertices of tri lie on the same size (positive/negative) of pl,
    // there is no intersection
    if (n_points_on_positive == 3 || n_points_on_negative == 3)
      return Object();

    // if all vertices of tri lie on pl then we return tri
    if (n_points_on_plane == 3)
       return make_object(tri);

    // if 2 vertices lie on pl, then return the segment between them
    if (n_points_on_plane == 2)
    {
      int point_idx1 = points_on_plane[0], point_idx2 = points_on_plane[1];
      return make_object (Segment_3(tri.vertex(point_idx1),
                                    tri.vertex(point_idx2)));
    }

    // if only 1 lie on pl, should check the segment opposite to it on tri
    if (n_points_on_plane == 1)
    {
      int point_on_plane_idx = points_on_plane[0];

      // if the other 2 vertices are on the same side of pl,
      // then the answer is just this vertex
      if (n_points_on_negative == 2 || n_points_on_positive == 2)
        return make_object(tri.vertex(point_on_plane_idx));

      // now it is known that one vertex is on pl, and the segment of tri
      // opposite to it should intersect pl

      // the segment of tri opposite of tri[point_on_plane_idx]
      Segment_3 tri_segment(tri.vertex(point_on_plane_idx+1),
                            tri.vertex(point_on_plane_idx+2));

      Object inter_result = k.intersect_3_object()(pl, tri_segment);
      Point_3 inter_point;
      CGAL_assertion( k.assign_3_object()(inter_point, inter_result) );
      k.assign_3_object()(inter_point, inter_result);

      // create the resulting segment
      // (between tri[point_on_plane_idx] and inter_point)
      return make_object(Segment_3(tri.vertex(point_on_plane_idx), 
                                   inter_point));

    }

    CGAL_assertion( n_points_on_plane == 0 );
    CGAL_assertion( n_points_on_positive + n_points_on_negative == 3 );
    CGAL_assertion( n_points_on_positive != 0 );
    CGAL_assertion( n_points_on_negative != 0 );

    // now it known that there is an intersection between 2 segments of tri
    // and pl, it is also known which segments are those.
    Point_3 inter_points[2];
    int pos_it, neg_it, n_inter_points = 0;
    for(pos_it = 0; pos_it < n_points_on_positive; ++pos_it)
      for(neg_it = 0; neg_it < n_points_on_negative; ++neg_it)
      {
        Segment_3 seg(tri.vertex(points_on_positive[pos_it]),
                      tri.vertex(points_on_negative[neg_it]));
        Object inter_result = k.intersect_3_object()(pl, seg);
        Point_3 inter_point;
        // the result of the intersection must be a point
        CGAL_assertion( k.assign_3_object()(inter_point, inter_result) );
        k.assign_3_object()(inter_point, inter_result);
        inter_points[n_inter_points++] = inter_point;
      }

    CGAL_assertion( n_inter_points == 2 );
    return make_object(Segment_3(inter_points[0], inter_points[1]));
  }

  // compare the value of s1 in p1 to the value of s2 in p2
  Comparison_result
  compare_z(const Point_2& p1,
            const Xy_monotone_surface_3& s1,
            const Point_2& p2,
            const Xy_monotone_surface_3& s2)
  {
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
  envelope_point_of_surface(const Point_2& p,
                            const Xy_monotone_surface_3& s) const
  {
    CGAL_precondition(s.is_xy_monotone());
    CGAL_precondition(is_defined_over_object()(p, s));

    Point_3 point(p.x(), p.y(), 0);

    // Compute the intersetion between the vertical line and the given surfaces
    if (s.is_segment())
      return envelope_point_of_segment(point, s);
    else
    {
      // s is a non-vertical triangle
      CGAL_assertion(!s.is_vertical());

      // Construct a vertical line passing through point
      Kernel k;
      Direction_3 dir (0, 0, 1);
      Line_3      vl = k.construct_line_3_object() (point, dir);

      const Plane_3& plane = s.plane();
      Object    res = k.intersect_3_object()(plane, vl);
      CGAL_assertion(!res.is_empty());
      Point_3 ip;
      CGAL_assertion(k.assign_3_object()(ip, res));
      k.assign_3_object()(ip, res);

      return ip;
    }
  }

  // find the envelope point of the surface over the given point
  // precondition: the surface is defined in point and is a segment
  Point_3 envelope_point_of_segment(const Point_3& point,
                                    const Xy_monotone_surface_3& s) const
  {
    Kernel k;
    CGAL_precondition(s.is_segment());
    CGAL_precondition(is_defined_over_object()(project(point), s));

    // this is the vertical plane through the segment
    const Plane_3& plane = s.plane();

    // Construct a vertical line passing through point
    Direction_3 dir (0, 0, 1);
    Line_3      vl = k.construct_line_3_object() (point, dir);
    // we need 2 points on this line, to be transformed to 2d,
    // and preserve the direction of the envelope
    Point_3 vl_point1 = k.construct_point_on_3_object()(vl, 0),
            vl_point2 = k.construct_point_on_3_object()(vl, 1);

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
    Object inter_obj = k.intersect_2_object()(seg, l);
    Point_2 inter_point;
    CGAL_assertion_code(bool is_inter_point =)
    k.assign_2_object()(inter_point, inter_obj);
    CGAL_assertion(is_inter_point);
    return plane.to_3d(inter_point);
  }

  Point_2 construct_middle_point(const Point_2& p1, const Point_2& p2) const
  {
    Kernel k;
    return k.construct_midpoint_2_object()(p1, p2);
  }

  Point_2 construct_middle_point(const X_monotone_curve_2& cv) const
  {
    Kernel k;
    return k.construct_midpoint_2_object()(cv.source(), cv.target());
  }

  /***************************************************************************/
  // for vertical decomposition
  /***************************************************************************/
  
  class Construct_vertical_2
  {
  public:
    X_monotone_curve_2 operator()(const Point_2& p1, const Point_2& p2) const
    {
      return X_monotone_curve_2(p1, p2);

    }
  };

  /*! Get a Construct_vertical_2 functor object. */
  Construct_vertical_2 construct_vertical_2_object() const
  {
    return Construct_vertical_2();
  }
 

  Point_2 vertical_ray_shoot_2 (const Point_2& pt,
                                const X_monotone_curve_2& cv) const
  {
    CGAL_precondition(!cv.is_vertical());

    typename Kernel::Segment_2 seg = cv;
    Kernel k;
    // If the curve contains pt, return it.
    if (k.has_on_2_object() (seg, pt))
      return (pt);

    // Construct a vertical line passing through pt.
    typename Kernel::Direction_2  dir (0, 1);
    typename Kernel::Line_2        vl = k.construct_line_2_object() (pt, dir);

    // Compute the intersetion between the vertical line and the given curve.
    Object    res = k.intersect_2_object()(seg, vl);
    Point_2   ip;
    bool      ray_shoot_successful = k.assign_2_object()(ip, res);

    if (! ray_shoot_successful)
      CGAL_assertion (ray_shoot_successful);

    return (ip);
  }
};


/*!
 * \class A representation of a triangle, as used by the 
 * Env_triangle_traits_3 traits-class.
 */
template <class Kernel_>
class Env_triangle_3 :
    public Env_triangle_traits_3<Kernel_>::_Triangle_cached_3
{
  typedef Kernel_                                                  Kernel;
  typedef typename Kernel::Triangle_3                              Triangle_3;
  typedef typename Kernel::Point_3                                 Point_3;
  typedef typename Kernel::Plane_3                                 Plane_3;
  typedef typename Kernel::Segment_3                               Segment_3;

  typedef typename Env_triangle_traits_3<Kernel>::_Triangle_cached_3
                                                                   Base;

public:

  /*!
   * Default constructor.
   */
  Env_triangle_3() :
    Base()
  {}

  /*!
   * Constructor from a "kernel" triangle.
   * \param seg The segment.
   */
  Env_triangle_3(const Triangle_3& tri) :
    Base(tri)
  {}

  /*!
   * Construct a triangle from 3 end-points.
   * \param p1 The first point.
   * \param p2 The second point.
   * \param p3 The third point.
   */
    Env_triangle_3(const Point_3 &p1, const Point_3 &p2, const Point_3 &p3) :
      Base(p1, p2, p3)
  {}

  /*!
   * Construct a triangle from a plane and 3 end-points.
   * \param pl The supporting plane.
   * \param p1 The first point.
   * \param p2 The second point.
   * \param p3 The third point.
   * \pre All points must be on the supporting plane.
   */
  Env_triangle_3(const Plane_3& pl,
                 const Point_3 &p1,
                 const Point_3 &p2,
                 const Point_3 &p3) :
    Base(pl, p1, p2, p3)

  {}

  /*!
   * Construct a segment from 2 end-points.
   * \param p1 The first point.
   * \param p2 The second point.
   */
  Env_triangle_3(const Point_3 &p1, const Point_3 &p2) :
    Base(p1, p2)
  {}

  /*!
   * Cast to a triangle.
   */
  operator Triangle_3() const
  {
    return (Triangle_3(this->vertex(0), this->vertex(1), this->vertex(2)));
  }

  /*!
   * Cast to a segment (only when possible).
   */
  operator Segment_3() const
  {
    CGAL_precondition(this->is_segment());
    return (Segment_3(this->vertex(0), this->vertex(1)));
  }

  /*!
   * Create a bounding box for the triangle.
   */
  Bbox_3 bbox() const
  {
    Triangle_3 tri(this->vertex(0), this->vertex(1), this->vertex(2));
    return (tri.bbox());
  }
};

template <class Kernel>
bool
operator<(const Env_triangle_3<Kernel> &a,
          const Env_triangle_3<Kernel> &b)
{
  if (a.vertex(0) < b.vertex(0))
    return true;
  if (a.vertex(0) > b.vertex(0))
    return false;
  if (a.vertex(1) < b.vertex(1))
    return true;
  if (a.vertex(1) > b.vertex(1))
    return false;
  if (a.vertex(2) < b.vertex(2))
    return true;
  if (a.vertex(2) > b.vertex(2))
    return false;

  return false;
}
template <class Kernel>
bool
operator==(const Env_triangle_3<Kernel> &a,
           const Env_triangle_3<Kernel> &b)
{
  return (a.vertex(0) == b.vertex(0) &&
          a.vertex(1) == b.vertex(1) &&
          a.vertex(2) == b.vertex(2));
}

/*!
 * Exporter for the triangle class used by the traits-class.
 */
template <class Kernel, class OutputStream>
OutputStream& operator<< (OutputStream& os, const Env_triangle_3<Kernel>& tri)
{
  os << static_cast<typename Kernel::Triangle_3>(tri);
  if (tri.is_segment())
    os << " (segment)";
  return (os);
}

/*!
 * Importer for the triangle class used by the traits-class.
 */
template <class Kernel, class InputStream>
InputStream& operator>> (InputStream& is, Env_triangle_3<Kernel>& tri)
{
  typename Kernel::Triangle_3   kernel_tri;
  is >> kernel_tri;
  tri = kernel_tri;
  return (is);
}

} //namespace CGAL

#endif 
