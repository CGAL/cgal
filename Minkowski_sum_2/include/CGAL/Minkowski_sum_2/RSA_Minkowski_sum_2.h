// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source: $
// $Revision$ $Date$
// $Name:  $
//
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_RSA_MINKOWSKI_SUM_2_H
#define CGAL_RSA_MINKOWSKI_SUM_2_H

#include <CGAL/Polygon_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Minkowski_sum_2/Arr_labeled_traits_2.h>
#include <CGAL/Minkowski_sum_2/Union_of_curve_cycles_2.h>
#include <list>

CGAL_BEGIN_NAMESPACE

/*! \class
 * A base class for computing the Minkowski sum of a linear polygons with
 * a polygon with circular arcs, which represents the rotational swept area
 * of some linear polygon.
 */
template <class ConicTraits_, class Container_>
class RSA_Minkowski_sum_2
{
private:
  
  // Linear kernel types:
  typedef ConicTraits_                                  Conic_traits_2;
  typedef typename Conic_traits_2::Rat_kernel           Rat_kernel;
  typedef typename Rat_kernel::FT                       Rational;
  typedef typename Rat_kernel::Point_2                  Rat_point_2;
  typedef typename Rat_kernel::Segment_2                Rat_segment_2;
  typedef typename Rat_kernel::Circle_2                 Rat_circle_2;
  typedef typename Rat_kernel::Vector_2                 Rat_vector_2;
  typedef typename Rat_kernel::Direction_2              Rat_direction_2;
  typedef Polygon_2<Rat_kernel, Container_>             Rat_polygon_2;

  // Linear kernel functors:
  typedef typename Rat_kernel::Equal_2                 Equal_2;
  typedef typename Rat_kernel::Compare_xy_2            Compare_xy_2;
  typedef typename Rat_kernel::Construct_center_2      Construct_center_2;
  typedef typename Rat_kernel::Construct_circle_2      Construct_circle_2;
  typedef typename Rat_kernel::Compute_squared_radius_2
                                                       Compute_sqr_radius_2;
  typedef typename Rat_kernel::Compute_squared_distance_2
                                                       Compute_sqr_distance_2;
  typedef typename Rat_kernel::Construct_vector_2      Construct_vector_2;
  typedef typename Rat_kernel::Construct_perpendicular_vector_2
                                                       Construct_perp_vector_2;
  typedef typename Rat_kernel::Construct_direction_2   Construct_direction_2;
  typedef typename Rat_kernel::Counterclockwise_in_between_2 Ccw_in_between_2;
  typedef typename Rat_kernel::Construct_translated_point_2  Translate_point_2;

  // Circle/segment traits types:
  typedef Gps_circle_segment_traits_2<Rat_kernel>        Gps_circ_traits_2;
  typedef typename Gps_circ_traits_2::Point_2            Circ_point_2;
  typedef typename Gps_circ_traits_2::X_monotone_curve_2 Circ_segment_2;
  typedef typename Gps_circ_traits_2::Polygon_2          Circ_polygon_2;

  // Conic traits types:
  typedef typename Conic_traits_2::Alg_kernel           Alg_kernel;
  typedef typename Conic_traits_2::Nt_traits            Nt_traits;
  typedef typename Alg_kernel::FT                       Algebraic;
  typedef typename Alg_kernel::Point_2                  Alg_point_2;
  typedef typename Conic_traits_2::Curve_2              Curve_2;
  typedef typename Conic_traits_2::X_monotone_curve_2   X_monotone_curve_2;

  typedef typename Alg_kernel::Compare_xy_2             Alg_compare_xy_2;

  typedef Gps_traits_2<Conic_traits_2>                  Gps_traits_2;
 
public:

  typedef typename Gps_traits_2::Polygon_2              Sum_polygon_2;
  
private:

  // Polygon-related types:
  typedef typename Rat_polygon_2::Vertex_const_circulator
                                                         Rat_vertex_circulator;
  typedef typename Circ_polygon_2::Curve_const_iterator  Circ_edge_iterator;

  // Labeled traits:
  typedef Arr_labeled_traits_2<Conic_traits_2>           Labeled_traits_2; 
  typedef typename Labeled_traits_2::X_monotone_curve_2  Labeled_curve_2;
  typedef std::list<Labeled_curve_2>                     Curves_list;

  typedef Union_of_curve_cycles_2<Labeled_traits_2,
                                  Sum_polygon_2>         Union_2;

  // Data members:
  Nt_traits               nt_traits;

  Equal_2                 f_equal;
  Construct_center_2      f_center;
  Construct_circle_2      f_circle;
  Compute_sqr_radius_2    f_sqr_radius;
  Compute_sqr_distance_2  f_sqr_distance;
  Construct_vector_2      f_vector;
  Construct_perp_vector_2 f_perp_vector;
  Construct_direction_2   f_direction;
  Ccw_in_between_2        f_ccw_in_between;
  Translate_point_2       f_add;
  Compare_xy_2            f_compare_xy;
  Alg_compare_xy_2        f_compare_xy_alg;

public:

  /*! Default constructor. */
  RSA_Minkowski_sum_2 ()
  {
    // Obtain the rational-kernel functors.
    Rat_kernel                ker;

    f_equal = ker.equal_2_object();
    f_center = ker.construct_center_2_object();
    f_circle = ker.construct_circle_2_object();
    f_sqr_radius = ker.compute_squared_radius_2_object();
    f_sqr_distance = ker.compute_squared_distance_2_object();
    f_vector = ker.construct_vector_2_object();
    f_perp_vector = ker.construct_perpendicular_vector_2_object();
    f_direction = ker.construct_direction_2_object();
    f_ccw_in_between = ker.counterclockwise_in_between_2_object();
    f_add = ker.construct_translated_point_2_object(); 
    f_compare_xy = ker.compare_xy_2_object();

    // Obtain the algebraic-kernel functors.
    Alg_kernel                alg_ker;

    f_compare_xy_alg = alg_ker.compare_xy_2_object();
  }

  /*!
   * Compute the Minkowski sum of a linear polygon and a polygon that contains
   * circular arcs, which represents the rotational swept-area of a linear
   * polygon.
   * Note that as the input polygons may not be convex, the Minkowski sum may
   * not be a simple polygon. The result is therefore represented as
   * the outer boundary of the Minkowski sum (which is always a simple polygon)
   * and a container of simple polygons, representing the holes inside this
   * polygon.
   * \param pgn1 The linear polygon.
   * \param pgn2 The polygon with circular arcs.
   * \param sum_bound Output: A polygon respresenting the outer boundary
   *                          of the Minkowski sum.
   * \param sum_holes Output: An output iterator for the holes in the sum,
   *                          represented as simple polygons.
   * \pre Both input polygons are simple.
   * \pre The value-type of the output iterator is Sum_polygon_2.
   * \return A past-the-end iterator for the holes in the sum.
   */
  template <class OutputIterator>
  OutputIterator operator() (const Rat_polygon_2& pgn1,
                             const Circ_polygon_2& pgn2,
                             Sum_polygon_2& sum_bound,
                             OutputIterator sum_holes) const
  {
    CGAL_precondition (pgn1.is_simple());

    // Compute the convolution of the two polygons.
    Curves_list      conv_cycle;

    _convolution_cycle (pgn1, pgn2,
                        1,           // Fictitious cycle ID.
                        conv_cycle);

    // Compute the union of the cycle(s) that represent the Minkowski sum.
    Union_2     unite;

    sum_holes = unite (conv_cycle.begin(), conv_cycle.end(),
                       sum_bound, sum_holes);

    return (sum_holes);
  }

protected:

  /*!
   * Compute the curves that constitute the convolution cycle(s) of a given
   * linear polygon with another polygon with circular arcs.
   * \param pgn1 The linear polygon.
   * \param pgn2 The polygon with circular arcs.
   * \param cycle_id The index of the cycle.
   * \param cycle Output: The curves that consitute the convolution cycle(s).
   */
  void _convolution_cycle (const Rat_polygon_2& pgn1,
                           const Circ_polygon_2& pgn2,
                           unsigned int cycle_id,
                           Curves_list& cycle) const
  {
    // Go over the vertices of pgn2 (at each iteration, the vertex we consider
    // is the target of prev2 and the source of curr2).
    const bool            forward1 = (pgn1.orientation() == 
                                      CGAL::COUNTERCLOCKWISE);
    Rat_vertex_circulator first1 = pgn1.vertices_circulator();
    Rat_vertex_circulator curr1, next1;
    Rat_direction_2       dir_curr1;
    Circ_edge_iterator    begin2 = pgn2.curves_begin();
    Circ_edge_iterator    end2 = pgn2.curves_end();
    Circ_edge_iterator    prev2, next2;
    Rat_point_2           p2;
    Rat_direction_2       dir_prev2, dir_next2;
    bool                  equal_dirs;
    bool                  shift_edge;
    unsigned int          xcv_index = 0;
    X_monotone_curve_2    shifted_xcv;
    bool                  dir_right;

    prev2 = end2; --prev2;
    next2 = begin2;
    while (next2 != end2)
    {
      // The current vertex in pgn2 must have rational coordinates.
      CGAL_assertion (prev2->target().x().is_rational());
      CGAL_assertion (prev2->target().y().is_rational());
      
      p2 = Rat_point_2 (prev2->target().x().alpha(),
                        prev2->target().y().alpha());

      // Get the directions of the edges around the current vertex in pgn2.
      dir_prev2 = _direction (*prev2, p2);
      dir_next2 = _direction (*next2, p2);

      equal_dirs = f_equal (dir_prev2, dir_next2);

      // Go over the edges of pgn1.
      next1 = curr1 = first1; 
      do
      {
        if (forward1)
          ++next1;
        else
          --next1;

        // Compute the direction of the current edge.
        dir_curr1 =  f_direction (f_vector (*curr1, *next1));

        if (! equal_dirs)
        {
          // Check if the current edge is between the two previously computed
          // directions.
          shift_edge = f_ccw_in_between (dir_curr1, dir_prev2, dir_next2);
        }
        else
        {
          // Check if the direction of the edge equals the two previously
          // computed directions, which are equal to one another.
          shift_edge = f_equal (dir_curr1, dir_prev2);
        }

        if (shift_edge)
        {
          // Shift the current edge in pgn1 by the current vertex in pgn2.
          shifted_xcv = _shift_segment (*curr1, *next1,
                                        p2,
                                        dir_right);

          cycle.push_back (Labeled_curve_2 (shifted_xcv,
                                            X_curve_label (dir_right,
                                                           cycle_id,
                                                           2*xcv_index,
                                                           false)));
          xcv_index++;
        }

        curr1 = next1;
      } while (curr1 != first1);

      // Move to the next vertex of pgn2.
      prev2 = next2;
      ++next2;
    }

    // Go over the vertices of pgn1.
    Rat_vertex_circulator prev1;
    Rat_direction_2       dir_prev1, dir_next1;
    Circ_edge_iterator    curr2;
    Rat_direction_2       dir_curr2;
    Rat_point_2           s2, t2;
    Rat_direction_2       dir_s2, dir_t2;
    bool                  is_between_s2, is_between_t2;
    bool                  is_between_prev1, is_between_next1;
    Alg_point_2           ps, pt;
    bool                  shift_arc, shift_prev1, shift_next1;

    prev1 = pgn1.vertices_circulator();
    if (forward1)
      --prev1;
    else
      ++prev1;

    next1 = curr1 = first1; 
    do
    {
      if (forward1)
        ++next1;
      else
        --next1;

      // Compute the direction of the two edges incident to the current vertex.
      dir_prev1 = f_direction (f_vector (*prev1, *curr1));
      dir_next1 = f_direction (f_vector (*curr1, *next1));

      // Go over all edges of pgn2.
      for (curr2 = begin2; curr2 != end2; ++curr2)
      {
        // Act according to the edge type.
        if (curr2->is_linear())
        {
          // Compute the direction of the linear edge.
          dir_curr2 = f_direction (curr2->supporting_line());
          
          // Check if the current edge is between the two previously computed
          // directions.
          if (f_ccw_in_between (dir_curr1, dir_prev2, dir_next2) ||
              f_equal (dir_curr1, dir_next2))
          {
            // Shift the current edge in pgn1 by the current vertex in pgn2.
            s2 = Rat_point_2 (curr2->source().x().alpha(),
                              curr2->source().y().alpha());
            
            t2 = Rat_point_2 (curr2->target().x().alpha(),
                              curr2->target().y().alpha());

            shifted_xcv = _shift_segment (s2, t2,
                                          *curr1,
                                          dir_right);

            cycle.push_back (Labeled_curve_2 (shifted_xcv,
                                              X_curve_label (dir_right,
                                                             cycle_id,
                                                             2*xcv_index,
                                                             false)));
            xcv_index++;
          }
        }
        else
        {
          // The current edge is a circular arc: compute the directions of
          // the arc tagents at its endpoints.
          s2 = Rat_point_2 (curr2->source().x().alpha(),
                            curr2->source().y().alpha());
          dir_s2 = _direction (*curr2, s2);

          t2 = Rat_point_2 (curr2->target().x().alpha(),
                            curr2->target().y().alpha());
          dir_t2 = _direction (*curr2, t2);

          CGAL_assertion (! f_equal (dir_s2, dir_t2));

          // Check if these two directions are between the directions of the
          // two edges incident to curr1.
          is_between_s2 = f_ccw_in_between (dir_s2, dir_prev1, dir_next1) ||
                          f_equal (dir_s2, dir_next1);
          
          is_between_t2 = f_ccw_in_between (dir_t2, dir_prev1, dir_next1) ||
                          f_equal (dir_t2, dir_next1);

          // Act according to the result.
          shift_arc = shift_prev1 = shift_next1 = false;

          if (is_between_s2 && is_between_t2)
          {
            // We should add the entire arc, shifted by curr1, to the cycle.
            ps = _convert_point (curr2->source());
            pt = _convert_point (curr2->target());
            shift_arc = true;
          }
          else if (is_between_s2)
          {
            // We should split the arc and shift the first portion. We should
            // also shift the segment after curr1 by the split point.
            ps = _convert_point (curr2->source());
            pt = _point_on_circle (curr2->supporting_circle(),
                                   *curr1, *next1);
            shift_arc = true;
            shift_next1 = true;
          }
          else if (is_between_t2)
          {
            // We should split the arc and shift the second portion. We should
            // also shift the segment before curr1 by the split point.
            ps = _point_on_circle (curr2->supporting_circle(),
                                   *prev1, *curr1);
            pt = _convert_point (curr2->target());

            shift_arc = true;
            shift_prev1 = true;
          }
          else
          {
            // Either the two directions of the edges around curr1 are both
            // between the two directions at the circular arc's endpoint,
            // or both are not (in the latter case, we do nothing).
            is_between_prev1 = f_ccw_in_between (dir_prev1, dir_s2, dir_t2) ||
                               f_equal (dir_prev1, dir_t2);

            is_between_next1 = f_ccw_in_between (dir_next1, dir_s2, dir_t2) ||
                               f_equal (dir_next1, dir_t2);

            if (is_between_prev1 && is_between_next1)
            {
              ps = _point_on_circle (curr2->supporting_circle(),
                                     *prev1, *curr1);
              pt = _point_on_circle (curr2->supporting_circle(),
                                     *curr1, *next1);
              
              shift_arc = true;
              shift_prev1 = true;
              shift_next1 = true;
            }
            else
            {
              CGAL_assertion (! is_between_prev1 && ! is_between_next1);
            }
          }

          // If necessary, shifted the circular arc by curr1.
          if (shift_arc)
          {
            shifted_xcv = _shift_circular_arc (curr2->supporting_circle(),
                                               ps, pt,
                                               *curr1,
                                               dir_right);

            cycle.push_back (Labeled_curve_2 (shifted_xcv,
                                              X_curve_label (dir_right,
                                                             cycle_id,
                                                             2*xcv_index,
                                                             false)));
            xcv_index++;
          }

          // If necessary, shift the polygon edge before curr1 by ps.
          if (shift_prev1)
          {
            shifted_xcv = _shift_segment (*prev1, *curr1,
                                          ps,
                                          dir_right);

            cycle.push_back (Labeled_curve_2 (shifted_xcv,
                                              X_curve_label (dir_right,
                                                             cycle_id,
                                                             2*xcv_index,
                                                             false)));
            xcv_index++;
          }

          // If necessary, shift the polygon edge after curr1 by pt.
          if (shift_next1)
          {
            shifted_xcv = _shift_segment (*curr1, *next1,
                                          pt,
                                          dir_right);

            cycle.push_back (Labeled_curve_2 (shifted_xcv,
                                              X_curve_label (dir_right,
                                                             cycle_id,
                                                             2*xcv_index,
                                                             false)));
            xcv_index++;
          }
        }
      }

      // Move to the next vertex of pgn1.
      prev1 = curr1;
      curr1 = next1;
      
    } while (curr1 != first1);

    return;
  }

  /*!
   * Compute the direction of a curve at a given point.
   * \param cv The curve (a line segment or a circular arc).
   * \param p The point (either its source or its target).
   * \return The direction.
   */
  Rat_direction_2 _direction (const Circ_segment_2& cv,
                              const Rat_point_2& p) const
  {
    // Act according to the curve type.
    if (cv.is_linear())
    {
      // The edge is a line segment, having a constant direction.
      return (f_direction (cv.supporting_line()));
    }
    else
    {
      // The curve is a circular arc. We compute the direction of the
      // tangent to its supporting circle at q, which is the direction of the
      // vector (q - c) rotated counterclockwise by 90 degrees, where c is the
      // center of the supporting circle.
      Rat_circle_2    circ = cv.supporting_circle();
      Rat_vector_2    vec = f_vector (f_center (circ), p);
      
      return (f_direction (f_perp_vector (vec, CGAL::COUNTERCLOCKWISE)));
    }
  }

  /*!
   * Convert a point with one-root coordinates to an algebraic point.
   */
  Alg_point_2 _convert_point (const Circ_point_2& p) const
  {
    Algebraic      x = nt_traits.convert (p.x().alpha());

    if (! p.x().is_rational())
    {
      x += nt_traits.convert (p.x().beta()) * 
           nt_traits.sqrt (nt_traits.convert (p.x().gamma()));
    }

    Algebraic      y = nt_traits.convert (p.y().alpha());

    if (! p.y().is_rational())
    {
      y += nt_traits.convert (p.y().beta()) * 
           nt_traits.sqrt (nt_traits.convert (p.y().gamma()));
    }

    return (Alg_point_2 (x, y));
  }

  /*!
   * Compute a point on the given circle, where the tangent to the circle
   * has the same direction as the given vector.
   * \param circ The circle.
   * \param ps The source point of the vector.
   * \param pt The target point of the vector.
   * \return The computed point.
   */
  Alg_point_2 _point_on_circle (const Rat_circle_2& circ,
                                const Rat_point_2& ps,
                                const Rat_point_2& pt) const
  {
    // The angle theta between the vector v and the x-axis is given by:
    //
    //                 delta_y                      delta_x
    //   sin(alpha) = ---------       cos(alpha) = ---------
    //                   len                          len
    //
    const Rat_point_2     pc = f_center (circ);
    const Rational        sqr_r = f_sqr_radius (circ);
    const Rational        delta_x = pt.x() - ps.x();
    const Rational        delta_y = pt.y() - ps.y();
    const Rational        sqr_len = f_sqr_distance (ps, pt);

    // To compute the desired point, we have to translate the circle
    // center by r in a direction that forms and angle of (alpha - PI/2)
    // with the x-axis, and we have:
    //
    //   trans_x = r*cos(alpha - PI/2) = r*sin(alpha)
    //   trans_y = r*sin(alpha - PI/2) = -r*cos(alpha)
    const Algebraic       r_over_len = nt_traits.sqrt (nt_traits.convert
                                                       (sqr_r / sqr_len));
    const Algebraic       trans_x = nt_traits.convert (delta_y) * r_over_len;
    const Algebraic       trans_y = nt_traits.convert (-delta_x) * r_over_len;

    return (Alg_point_2 (nt_traits.convert (pc.x()) + trans_x, 
                         nt_traits.convert (pc.y()) + trans_y));
  }

  /*!
   * Shift the given segment by a point with rational coordinates.
   * \param ps The source point of the segment.
   * \param pt The target point of the segment.
   * \param q The shift point.
   * \param dir_right Output: Is the segment directed right.
   * \return The shifted segment, represented as an x-monotone conic curve.
   */
  X_monotone_curve_2 _shift_segment (const Rat_point_2& ps,
                                     const Rat_point_2& pt,
                                     const Rat_point_2& q,
                                     bool& dir_right) const
  {
    // Determine if the segment is directed left or right.
    Comparison_result  res = f_compare_xy (ps, pt);

    CGAL_assertion (res != CGAL::EQUAL);
    dir_right = (res == CGAL::SMALLER);

    // Construct the shifted segment.
    Rat_point_2  shift_ps = f_add (ps, f_vector(CGAL::ORIGIN, q));
    Rat_point_2  shift_pt = f_add (pt, f_vector(CGAL::ORIGIN, q));
    Curve_2      seg (Rat_segment_2 (shift_ps, shift_pt));

    return (seg);
  }

  /*!
   * Shift the given segment by a point with algebraic coordinates.
   * \param ps The source point of the segment.
   * \param pt The target point of the segment.
   * \param q The shift point.
   * \param dir_right Output: Is the segment directed right.
   * \return The shifted segment, represented as an x-monotone conic curve.
   */
  X_monotone_curve_2 _shift_segment (const Rat_point_2& ps,
                                     const Rat_point_2& pt,
                                     const Alg_point_2& q,
                                     bool& dir_right) const
  {
    // Determine if the segment is directed left or right.
    Comparison_result  res = f_compare_xy (ps, pt);

    CGAL_assertion (res != CGAL::EQUAL);
    dir_right = (res == CGAL::SMALLER);

    // Construct the shifted segment.
    Alg_point_2  shift_ps = Alg_point_2 (nt_traits.convert (ps.x()) + q.x(), 
                                         nt_traits.convert (ps.y()) + q.y());
    Alg_point_2  shift_pt = Alg_point_2 (nt_traits.convert (pt.x()) + q.x(), 
                                         nt_traits.convert (pt.y()) + q.y());

    // The supprting line a*x + b*y + c = 0 connecting the two shifted points
    // is given by:
    const Algebraic  a = nt_traits.convert (ps.y() - pt.y());
    const Algebraic  b = nt_traits.convert (pt.x() - ps.x());
    const Algebraic  c = nt_traits.convert (ps.x()*pt.y() - ps.y()*pt.x()) -
                         (a*q.x() + b*q.y());
    
    return (X_monotone_curve_2 (a, b, c,
                                shift_ps, shift_pt));
  }

  /*!
   * Shift a circular arc, given by a rational supporting circle and two
   * endpoints, by a point with rational coordinates.
   * \param circ The supporting circle.
   * \param ps The source point of the arc.
   * \param pt The target point of the arc.
   * \param q The shift point.
   * \param dir_right Output: Is the segment directed right.
   * \return The shifted arc.
   */
  X_monotone_curve_2 _shift_circular_arc (const Rat_circle_2& circ,
                                          const Alg_point_2& ps,
                                          const Alg_point_2& pt,
                                          const Rat_point_2& q,
                                          bool& dir_right) const
  {
    // Determine if the segment is directed left or right.
    Comparison_result  res = f_compare_xy_alg (ps, pt);

    CGAL_assertion (res != CGAL::EQUAL);
    dir_right = (res == CGAL::SMALLER);

    // Shift the supporting circle.
    Rat_point_2  shift_pc = f_add (f_center (circ), f_vector(CGAL::ORIGIN, q));
    Rat_circle_2 shift_circ = f_circle (shift_pc,
                                        f_sqr_radius (circ));

    // Shift the two endpoint and construct the circular arc.
    Alg_point_2  shift_ps = Alg_point_2 (ps.x() + nt_traits.convert (q.x()), 
                                         ps.y() + nt_traits.convert (q.y()));
    Alg_point_2  shift_pt = Alg_point_2 (pt.x() + nt_traits.convert (q.x()), 
                                         pt.y() + nt_traits.convert (q.y()));
    Curve_2      arc (shift_circ, 
                      CGAL::COUNTERCLOCKWISE,
                      shift_ps, shift_pt);

    return (arc);
  }

};


CGAL_END_NAMESPACE

#endif
