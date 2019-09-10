// Copyright (c) 2006-2008  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein       <wein_r@yahoo.com>
//                 Andreas Fabri  <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>
//                 Efi Fogel      <efif@post.tau.ac.il>

#ifndef CGAL_APPROXIMATED_OFFSET_BASE_H
#define CGAL_APPROXIMATED_OFFSET_BASE_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Minkowski_sum_2/Labels.h>
#include <CGAL/Minkowski_sum_2/Arr_labeled_traits_2.h>
#include <CGAL/use.h>

namespace CGAL {

/*! \class
 * A base class for approximating the offset of a given polygon by a given
 * radius.
 */
template <class Kernel_, class Container_>
class Approx_offset_base_2
{
private:

  typedef Kernel_                                        Kernel;
  typedef typename Kernel::FT                            NT;

protected:

  typedef Kernel                                         Basic_kernel;
  typedef NT                                             Basic_NT;
  typedef CGAL::Polygon_2<Kernel, Container_>            Polygon_2;
  typedef CGAL::Polygon_with_holes_2<Kernel, Container_> Polygon_with_holes_2;

private:

  // Kernel types:
  typedef typename Kernel::Point_2                       Point_2;
  typedef typename Kernel::Line_2                        Line_2;

  // Polygon-related types:
  typedef typename Polygon_2::Vertex_circulator          Vertex_circulator;

  // Traits-class types:
  typedef Gps_circle_segment_traits_2<Kernel>            Traits_2;
  typedef typename Traits_2::CoordNT                     CoordNT;
  typedef typename Traits_2::Point_2                     Tr_point_2;
  typedef typename Traits_2::Curve_2                     Curve_2;
  typedef typename Traits_2::X_monotone_curve_2          X_monotone_curve_2;

protected:

  typedef typename Traits_2::Polygon_2                   Offset_polygon_2;

  typedef Arr_labeled_traits_2<Traits_2>                 Labeled_traits_2;
  typedef typename Labeled_traits_2::X_monotone_curve_2  Labeled_curve_2;

  // Data members:
  double        _eps;            // An upper bound on the approximation error.
  int           _inv_sqrt_eps;   // The inverse squared root of _eps.

public:

  /*!
   * Constructor.
   * \param eps An upper bound on the approximation error.
   */
  Approx_offset_base_2 (const double& eps) :
    _eps (eps)
  {
    CGAL_precondition (CGAL::sign (eps) == POSITIVE);

    _inv_sqrt_eps = static_cast<int> (1.0 / CGAL::sqrt (_eps));
    if (_inv_sqrt_eps <= 0)
      _inv_sqrt_eps = 1;
  }

protected:

  /*!
   * Compute curves that constitute the offset of a simple polygon by a given
   * radius, with a given approximation error.
   * \param pgn The polygon.
   * \param orient The orientation to traverse the vertices.
   * \param r The offset radius.
   * \param cycle_id The index of the cycle.
   * \param oi An output iterator for the offset curves.
   * \pre The value type of the output iterator is Labeled_curve_2.
   * \return A past-the-end iterator for the holes container.
   */
  template <class OutputIterator>
  OutputIterator _offset_polygon (const Polygon_2& pgn,
                                  CGAL::Orientation orient,
                                  const Basic_NT& r,
                                  unsigned int cycle_id,
                                  OutputIterator oi) const
  {
    // Prepare circulators over the polygon vertices.
    const bool            forward = (pgn.orientation() == orient);
    Vertex_circulator     first, curr, next, prev;

    first = pgn.vertices_circulator();
    curr = first;
    next = first;
    prev = first;

    if (forward)
      --prev;
    else
      ++prev;

    // Traverse the polygon vertices and edges and approximate the arcs that
    // constitute the single convolution cycle.
    NT              x1, y1;              // The source of the current edge.
    NT              x2, y2;              // The target of the current edge.
    NT              delta_x, delta_y;    // (x2 - x1) and (y2 - y1), resp.
    NT              abs_delta_x;
    NT              abs_delta_y;
    CGAL::Sign      sign_delta_x;        // The sign of (x2 - x1).
    CGAL::Sign      sign_delta_y;        // The sign of (y2 - y1).
    NT              sqr_d;               // The squared length of the edge.
    NT              err_bound;           // An approximation bound for d.
    NT              app_d;               // The apporximated edge length.
    NT              app_err;             // The approximation error.
    CGAL::Sign      sign_app_err;        // Its sign.
    NT              lower_tan_half_phi;
    NT              upper_tan_half_phi;
    NT              sqr_tan_half_phi;
    NT              sin_phi, cos_phi;
    Point_2         op1;                 // The approximated offset point
                                         // corresponding to (x1, y1).
    Point_2         op2;                 // The approximated offset point
                                         // corresponding to (x2, y2).
    Line_2          l1, l2;              // Lines tangent at op1 and op2.
    Object          obj;
    bool            assign_success;
    Point_2         mid_p;               // The intersection of l1 and l2.
    Point_2         first_op;            // op1 for the first edge visited.
    Point_2         prev_op;             // op2 for the previous edge.

    unsigned int                    curve_index = 0;
    X_monotone_curve_2              seg1, seg2;
    bool                            dir_right1 = false, dir_right2 = false;
    X_monotone_curve_2              seg_short;
    bool                            dir_right_short;
    int                             n_segments;

    Kernel                            ker;
    typename Kernel::Intersect_2      f_intersect = ker.intersect_2_object();
    typename Kernel::Construct_line_2 f_line = ker.construct_line_2_object();
    typename Kernel::Construct_perpendicular_line_2
                     f_perp_line = ker.construct_perpendicular_line_2_object();
    typename Kernel::Compare_xy_2     f_comp_xy = ker.compare_xy_2_object();
    typename Kernel::Orientation_2    f_orient = ker.orientation_2_object();

    Traits_2                        traits;
    std::list<Object>               xobjs;
    std::list<Object>::iterator     xobj_it;
    typename Traits_2::Make_x_monotone_2
                       f_make_x_monotone = traits.make_x_monotone_2_object();
    Curve_2                         arc;
    X_monotone_curve_2              xarc;

    do
    {
      // Get a circulator for the next vertex (in the proper orientation).
      if (forward)
        ++next;
      else
        --next;

      // Compute the vector v = (delta_x, delta_y) of the current edge,
      // and compute the squared edge length.
      x1 = curr->x();
      y1 = curr->y();
      x2 = next->x();
      y2 = next->y();

      delta_x = x2 - x1;
      delta_y = y2 - y1;
      sqr_d = CGAL::square (delta_x) + CGAL::square (delta_y);

      sign_delta_x = CGAL::sign (delta_x);
      sign_delta_y = CGAL::sign (delta_y);

      if (sign_delta_x == CGAL::ZERO)
      {
        CGAL_assertion (sign_delta_y != CGAL::ZERO);

        // The edge [(x1, y1) -> (x2, y2)] is vertical. The offset edge lies
        // at a distance r to the right if y2 > y1, and to the left if y2 < y1.
        if (sign_delta_y == CGAL::POSITIVE)
        {
          op1 = Point_2 (x1 + r, y1);
          op2 = Point_2 (x2 + r, y2);
        }
        else
        {
          op1 = Point_2 (x1 - r, y1);
          op2 = Point_2 (x2 - r, y2);
        }

        // Create the offset segment [op1 -> op2].
        seg1 = X_monotone_curve_2 (op1, op2);
        dir_right1 = (sign_delta_y == CGAL::POSITIVE);

        n_segments = 1;
      }
      else if (sign_delta_y == CGAL::ZERO)
      {
        // The edge [(x1, y1) -> (x2, y2)] is horizontal. The offset edge lies
        // at a distance r to the bottom if x2 > x1, and to the top if x2 < x1.
        if (sign_delta_x == CGAL::POSITIVE)
        {
          op1 = Point_2 (x1, y1 - r);
          op2 = Point_2 (x2, y2 - r);
        }
        else
        {
          op1 = Point_2 (x1, y1 + r);
          op2 = Point_2 (x2, y2 + r);
        }

        // Create the offset segment [op1 -> op2].
        seg1 = X_monotone_curve_2 (op1, op2);
        dir_right1 = (sign_delta_x == CGAL::POSITIVE);

        n_segments = 1;
      }
      else
      {
        abs_delta_x = (sign_delta_x == POSITIVE) ? delta_x : -delta_x;
        abs_delta_y = (sign_delta_y == POSITIVE) ? delta_y : -delta_y;

        // In this general case, the length d of the current edge is usually
        // an irrational number.
        // Compute the upper bound for the approximation error for d.
        // This bound is given by:
        //
        //                            d - |delta_y|
        //     bound = 2 * d * eps * ---------------
        //                              |delta_x|
        //
        // As we use floating-point arithmetic, if |delta_x| is small, then
        // it might be that to_double(|delta_y|) == to_double(d), hence we
        // have a 0 tolerance in the approximation bound. Luckily, because
        // of symmetry, we can rotate the scene by pi/2, and swap roles of
        // x and y. In fact, we do that in order to get a larger approximation
        // bound if possible.
        const double    dd = CGAL::sqrt (CGAL::to_double (sqr_d));
        const double    dabs_dx = CGAL::to_double (abs_delta_x);
        const double    dabs_dy = CGAL::to_double (abs_delta_y);
        double          derr_bound;

        if (dabs_dy < dabs_dx)
        {
          derr_bound = 2 * dd * _eps * (dd - dabs_dy) / dabs_dx;
        }
        else
        {
          derr_bound = 2 * dd * _eps * (dd - dabs_dx) / dabs_dy;
        }

        CGAL_assertion (derr_bound > 0);
        err_bound = NT (derr_bound);

        // Compute an approximation for d (the squared root of sqr_d).
        int             numer;
        int             denom = _inv_sqrt_eps;
        const int       max_int = (1 << (8*sizeof(int) - 2));

        numer = static_cast<int> (dd * denom + 0.5);
        if (numer > 0)
        {
            while (static_cast<double>(max_int) / denom < dd &&
                   numer > 0)
            {
              denom >>= 1;
              numer = static_cast<int> (dd * denom + 0.5);
            }
        }
        else if (numer == 0)
        {
            while (numer == 0)
            {
              denom <<= 1;
              if (denom > 0)
              {
                numer = static_cast<int> (dd * denom + 0.5);
              }
              else
              {
		// In case of overflow of denom
                numer = 1;
                denom = max_int;
              }
            }
        }
        else {// if numer < 0 (overflow)
          numer = max_int;
          denom = 1;
        }


        app_d = NT (numer) / NT (denom);
        app_err = sqr_d - CGAL::square (app_d);

        while (CGAL::compare (CGAL::abs (app_err),
                              err_bound) == CGAL::LARGER ||
               CGAL::compare (app_d, abs_delta_x) != LARGER ||
               CGAL::compare (app_d, abs_delta_y) != LARGER)
        {
          app_d = (app_d + sqr_d/app_d) / 2;
          app_err = sqr_d - CGAL::square (app_d);
        }

        sign_app_err = CGAL::sign (app_err);

        if (sign_app_err == CGAL::ZERO)
        {
          // In this case d is a rational number, and we should shift the
          // both edge endpoints by (r * delta_y / d, -r * delta_x / d) to
          // obtain the offset points op1 and op2.
          const NT   trans_x = r * delta_y / app_d;
          const NT   trans_y = r * (-delta_x) / app_d;

          op1 = Point_2 (x1 + trans_x, y1 + trans_y);
          op2 = Point_2 (x2 + trans_x, y2 + trans_y);

          seg1 = X_monotone_curve_2 (op1, op2);
          dir_right1 = (sign_delta_x == CGAL::POSITIVE);

          n_segments = 1;
        }
        else
        {
          // In case |x2 - x1| < |y2 - y1| (and phi is small) it is possible
          // that the approximation t' of t = tan(phi/2) is of opposite sign.
          // To avoid this problem, we symbolically rotate the scene by pi/2,
          // swapping roles between x and y. Thus, t is not close to zero, and
          // we are guaranteed to have: phi- < phi < phi+ .
          bool rotate_pi2 = false;

          if (CGAL::compare (CGAL::abs(delta_x),
                             CGAL::abs(delta_y)) == SMALLER)
          {
            rotate_pi2 = true;

            // We use the rotation matrix by pi/2:
            //
            //   +-       -+
            //   |  0  -1  |
            //   |  1   0  |
            //   +-       -+
            //
            // Thus, the point (x, y) is converted to (-y, x):
            NT tmp = x1;

            x1 = -y1;
            y1 = tmp;

            tmp = x2;
            x2 = -y2;
            y2 = tmp;

            // Swap the delta_x and delta_y values.
            tmp = delta_x;
            delta_x = -delta_y;
            delta_y = tmp;

            CGAL::Sign  tmp_sign = sign_delta_x;
            sign_delta_x = CGAL::opposite (sign_delta_y);
            sign_delta_y = tmp_sign;
          }

          // Act according to the sign of delta_x.
          if (sign_delta_x == CGAL::NEGATIVE)
          {
            // x1 > x2, so we take a lower approximation for the squared root.
            if (sign_app_err == CGAL::NEGATIVE)
              app_d = sqr_d / app_d;
          }
          else
          {
            // x1 < x2, so we take an upper approximation for the squared root.
            if (sign_app_err == CGAL::POSITIVE)
              app_d = sqr_d / app_d;
          }

          // If theta is the angle that the vector (delta_x, delta_y) forms
          // with the x-axis, the perpendicular vector forms an angle of
          // phi = theta - PI/2, and we can approximate tan(phi/2) from below
          // and from above using:
          lower_tan_half_phi = (app_d - delta_y) / (-delta_x);
          upper_tan_half_phi = (-delta_x) / (app_d + delta_y);

          // Translate (x1, y1) by (r*cos(phi-), r*sin(phi-)) and create the
          // first offset point.
          // If tan(phi/2) = t is rational, then sin(phi) = 2t/(1 + t^2)
          // and cos(phi) = (1 - t^2)/(1 + t^2) are also rational.
          sqr_tan_half_phi = CGAL::square (lower_tan_half_phi);
          sin_phi = 2 * lower_tan_half_phi / (1 + sqr_tan_half_phi);
          cos_phi = (1 - sqr_tan_half_phi) / (1 + sqr_tan_half_phi);

          if (! rotate_pi2)
          {
            op1 = Point_2 (x1 + r*cos_phi, y1 + r*sin_phi);
          }
          else
          {
            // In case of symbolic rotation by pi/2, we have to rotate the
            // translated point by -(pi/2), transforming (x, y) to (y, -x).
            op1 = Point_2 (y1 + r*sin_phi, -(x1 + r*cos_phi));
          }

          // Translate (x2, y2) by (r*cos(phi+), r*sin(phi+)) and create the
          // second offset point.
          sqr_tan_half_phi = CGAL::square (upper_tan_half_phi);
          sin_phi = 2 * upper_tan_half_phi / (1 + sqr_tan_half_phi);
          cos_phi = (1 - sqr_tan_half_phi) / (1 + sqr_tan_half_phi);

          if (! rotate_pi2)
          {
            op2 = Point_2 (x2 + r*cos_phi, y2 + r*sin_phi);
          }
          else
          {
            // In case of symbolic rotation by pi/2, we have to rotate the
            // translated point by -(pi/2), transforming (x, y) to (y, -x).
            op2 = Point_2 (y2 + r*sin_phi, -(x2 + r*cos_phi));
          }

          // Compute the line l1 tangent to the circle centered at (x1, y1)
          // with radius r at the approximated point op1.
          l1 = f_perp_line (f_line (*curr, op1), op1);

          // Compute the line l2 tangent to the circle centered at (x2, y2)
          // with radius r at the approximated point op2.
          l2 = f_perp_line (f_line (*next, op2), op2);

          // Intersect the two lines. The intersection point serves as a common
          // end point for the two line segments we are about to introduce.
          obj = f_intersect (l1, l2);

          assign_success = CGAL::assign (mid_p, obj);
          CGAL_assertion (assign_success);
          CGAL_USE(assign_success);

          // Andreas's assertions:
          CGAL_assertion( right_turn(*curr, *next, op2) );
          CGAL_assertion( angle(*curr, *next, op2) != ACUTE);
          CGAL_assertion( angle(op1, *curr, *next) != ACUTE);
          CGAL_assertion( right_turn(op1, *curr, *next) );

          // Create the two segments [op1 -> p_mid] and [p_min -> op2].
          seg1 = X_monotone_curve_2 (op1, mid_p);
          dir_right1 = (f_comp_xy (op1, mid_p) == CGAL::SMALLER);

          seg2 = X_monotone_curve_2 (mid_p, op2);
          dir_right2 = (f_comp_xy (mid_p, op2) == CGAL::SMALLER);

          n_segments = 2;
        }
      }

      if (curr == first) {
        // This is the first edge we visit -- store op1 for future use.
        first_op = op1;
      }
      else {
        CGAL::Orientation orient = f_orient (*prev, *curr, *next);
        if (orient == CGAL::COLLINEAR) {
          /* If the orientation is collinear, figure out whether it's a 180
           * turn. If so, assume that it is an antena that generates
           * a positive spike, and treat it as a left turn.
           * A complete solution would need to distinguish between positive
           * and negative spikes, and treat them as left and right turns,
           * respectively.
           */
          typename Kernel::Compare_x_2 f_comp_x = ker.compare_x_2_object();
          Comparison_result res1, res2;
          res1 = f_comp_x(*prev, *curr);
          if (res1 != CGAL::EQUAL)
            res2 = f_comp_x(*curr, *next);
          else {
            typename Kernel::Compare_y_2 f_comp_y = ker.compare_y_2_object();
            res1 = f_comp_y(*prev, *curr);
            res2 = f_comp_y(*curr, *next);
          }
          if (res1 != res2) orient = CGAL::LEFT_TURN;
        }

        // Connect the offset target point of the previous edge to the
        // offset source of the current edge.
        if (orient == CGAL::LEFT_TURN) {
          // Connect prev_op and op1 with a circular arc, whose supporting
          // circle is (x1, x2) with radius r.
          arc = Curve_2 (*curr, r, CGAL::COUNTERCLOCKWISE,
                         Tr_point_2 (prev_op.x(), prev_op.y()),
                         Tr_point_2 (op1.x(), op1.y()));

          // Subdivide the arc into x-monotone subarcs and insert them into the
          // convolution cycle.
          xobjs.clear();
          f_make_x_monotone (arc, std::back_inserter(xobjs));
          for (xobj_it = xobjs.begin(); xobj_it != xobjs.end(); ++xobj_it) {
            assign_success = CGAL::assign (xarc, *xobj_it);
            CGAL_assertion (assign_success);
            CGAL_USE(assign_success);
            *oi++ = Labeled_curve_2 (xarc,
                                     X_curve_label (xarc.is_directed_right(),
                                                    cycle_id, curve_index++));
          }
        }
        else if (orient == CGAL::RIGHT_TURN) {
          // In case the current angle between the previous and the current
          // edge is larger than pi/2, it not necessary to connect prev_op
          // and op1 by a circular arc (as the case above): it is sufficient
          // to shortcut the circular arc using a segment, whose sole purpose
          // is to guarantee the continuity of the convolution cycle (we know
          // this segment will not be part of the output offset or inset).
          seg_short = X_monotone_curve_2(prev_op, op1);

          dir_right_short = (f_comp_xy (prev_op, op1) == CGAL::SMALLER);
          *oi++ = Labeled_curve_2 (seg_short,
                                   X_curve_label (dir_right_short,
                                                  cycle_id, curve_index++));
        }
      }

      // Append the offset segment(s) to the convolution cycle.
      CGAL_assertion (n_segments == 1 || n_segments == 2);
      *oi++ = Labeled_curve_2 (seg1, X_curve_label (dir_right1,
                                                    cycle_id, curve_index++));

      if (n_segments == 2)
      {
        *oi++ = Labeled_curve_2 (seg2, X_curve_label (dir_right2,
                                                      cycle_id, curve_index++));
      }

      // Proceed to the next polygon vertex.
      prev_op = op2;
      prev = curr;
      curr = next;

    } while (curr != first);

    // Close the convolution cycle by creating the final circular arc,
    // centered at the first vertex.
    arc = Curve_2 (*first, r, CGAL::COUNTERCLOCKWISE,
                   Tr_point_2 (op2.x(), op2.y()),
                   Tr_point_2 (first_op.x(), first_op.y()));

    // Subdivide the arc into x-monotone subarcs and insert them to the
    // convolution cycle.
    xobjs.clear();
    f_make_x_monotone (arc, std::back_inserter(xobjs));

    xobj_it = xobjs.begin();
    while (xobj_it != xobjs.end())
    {
      assign_success = CGAL::assign (xarc, *xobj_it);
      CGAL_assertion (assign_success);
      CGAL_USE(assign_success);

      ++xobj_it;
      bool is_last = (xobj_it == xobjs.end());

      *oi++ = Labeled_curve_2 (xarc,
                               X_curve_label (xarc.is_directed_right(),
                                              cycle_id, curve_index++, is_last));
    }

    return (oi);
  }

};

} //namespace CGAL

#endif
