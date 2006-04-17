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

CGAL_BEGIN_NAMESPACE

/*! \class
 * A base class for computing the Minkowski sum of a linear polygons with
 * a polygon with circular arcs, which represents the rotational swept area
 * of some linear polygon.
 */
template <class ConicTraits_, class Container_, class DecompStrategy_>
class RSA_Minkowski_sum_2
{
private:
  
  // Linear kernel types:
  typedef ConicTraits_                                  Conic_traits_2;
  typedef typename Conic_traits_2::Rat_kernel           Rat_kernel;
  typedef typename Rat_kernel::FT                       Rational;
  typedef typename Rat_kernel::Point_2                  Rat_point_2;
  typedef typename Rat_kernel::Segment_2                Rat_segment_2;
  typedef typename Rat_kernel::Vector_2                 Rat_vector_2;
  typedef typename Rat_kernel::Direction_2              Rat_direction_2;
  typedef Polygon_2<Rat_kernel, Container>              Rat_polygon_2;

  // Linear kernel functors:
  typedef typename Rat_kernel::Equal_2                 Equal_2;
  typedef typename Rat_kernel::Compare_xy_2            Compare_xy_2;
  typedef typename Rat_kernel::Construct_center_2      Construct_center_2;
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

  typedef Gps_traits_2<Conic_traits_2>                  Gps_traits_2;
 
public:

  typedef typename Gps_traits_2::Polygon_2              Sum_polygon_2;
  
private:

  // Polygon-related types:
  typedef typename Polygon_2::Vertex_const_circulator    Rat_vertex_circulator;
  typedef typename Circ_polygon_2::Curve_const_iterator  Circ_edge_iterator;

  // Labeled traits:
  typedef Arr_labeled_traits_2<Conic_traits_2>           Labeled_traits_2; 
  typedef typename Labeled_traits_2::X_monotone_curve_2  Labeled_curve_2;

  // Data members:
  Nt_traits               nt_traits;

  Equal_2                 f_equal;
  Construct_center_2      f_center;
  Compute_sqr_radius_2    f_sqr_radius;
  Compute_sqr_distance_2  f_sqr_distance;
  Construct_vector_2      f_vector;
  Construct_perp_vector_2 f_perp_vector;
  Construct_direction_2   f_direction;
  Ccw_in_between_2        f_ccw_in_between;
  Translate_point_2       f_add;
  Compare_xy_2            f_compare_xy;

public:

  /*! Default constructor. */
  RSA_Minkowski_sum_2 ()
  {
    // Obtain the rational-kernel functors.
    Rat_kernel                ker;

    f_equal = ker.equal_2_object();
    f_center = ker.construct_center_2_object();
    f_sqr_radius = ker.compute_squared_radius_2_object();
    f_sqr_distance = ker.compute_squared_distance_2_object();
    f_vector = ker.construct_vector_2_object();
    f_perp_vector = ker.construct_perpendicular_vector_2_object();
    f_direction = ker.construct_direction_2_object();
    f_ccw_in_between = ker.counterclockwise_in_between_2_object();
    f_add = ker.construct_translated_point_2_object(); 
    f_compare_xy = ker.compare_xy_2_object();
  }

protected:

  /*!
   * Compute the curves that constitute the Minkowski sum of a convex linear
   * polygon with another polygon with circular arcs.
   * \param pgn1 The linear polygon.
   * \param pgn2 The polygon with circular arcs.
   * \param cycle_id The index of the cycle.
   * \param oi An output iterator for the offset curves.
   * \pre The value type of the output iterator is Labeled_curve_2.
   * \return A past-the-end iterator for the holes container.
   */
  template <class OutputIterator>
  OutputIterator _sum_with_convex (const Rat_polygon_2& pgn1,
                                   const Circ_polygon_2& pgn2,
                                   unsigned int cycle_id,
                                   OutputIterator oi) const
  {
    // Go over the vertices of pgn2 (at each iteration, the vertex we consider
    // is the target of prev2 and the source of curr2).
    const bool            forward1 = (pgn1.orientation() == COUNTERCLOCKWISE);
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
    Rat_point_2           ps, pt;
    unsigned int          xcv_index = 0;

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
          ps = f_add (*curr1, f_vector(CGAL::ORIGIN, p2));
          pt = f_add (*next1, f_vector(CGAL::ORIGIN, p2));
          
          res = f_compare_xy (ps, pt);
          CGAL_assertion (res != EQUAL);

          Curve_2    seg (Rat_segment_2 (ps, pt));
          cycle.push_back (Labeled_curve_2 (X_monotone_curve_2 (seg),
                                            X_curve_label ((res == SMALLER),
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
            ps = f_add (s2, f_vector(CGAL::ORIGIN, *curr1));
            
            t2 = Rat_point_2 (curr2->target.x().alpha(),
                              curr2->target().y().alpha());
            pt = f_add (t2, f_vector(CGAL::ORIGIN, *curr1));

            res = f_compare_xy (ps, pt);
            CGAL_assertion (res != EQUAL);

            Curve_2    seg (Rat_segment_2 (ps, pt));
            cycle.push_back (Labeled_curve_2 (X_monotone_curve_2 (seg),
                                              X_curve_label ((res == SMALLER),
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

          t2 = Rat_point_2 (curr2->target.x().alpha(),
                            curr2->target().y().alpha());
          dir_t2 = _direction (*curr2, s2);

          CGAL_assertion (! f_equal (dir_s2, dir_t2));

          // Check if these two directions are between the directions of the
          // two edges incident to curr1.
          is_between_s2 = f_ccw_in_between (dir_s2, dir_prev1, dir_next1) ||
                          f_equal (dir_s2, dir_next1);
          
          is_between_t2 = f_ccw_in_between (dir_t2, dir_prev1, dir_next1) ||
                          f_equal (dir_t2, dir_next1);

          // Act according to the result.
          if (is_between_s2 && is_between_t2)
          {
            // RWRW: (case (b))
          }
          else if (is_between_s2)
          {
            // RWRW: (case (c))
          }
          else if (is_between_t2)
          {
            // RWRW: (case (d)) 
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
              // RWRW: (case (e))
            }
            else
            {
              CGAL_assertion (! is_between_prev1 && ! is_between_next1);
            }
          }
        }
      }

      // Move to the next vertex of pgn1.
      prev1 = curr1;
      curr1 = next1;
      
    } while (curr1 != first1);

    return (oi);
  }

  /*!
   * Compute the direction of a curve at a given point.
   * \param cv The curve (a line segment or a circular arc).
   * \param p The point (either its source or its target).
   * \return The direction.
   */
  Direction_2 _direction (const Circ_segment_2& cv,
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
      
      return (f_direction (f_perp_vec (vec, COUNTERCLOCKWISE)));
    }
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
   * \return The shifted segemnt, represented as an x-monotone conic curve.
   */
  X_monotone_curve_2 _shift_segment (const Rat_point_2& ps,
                                     const Rat_point_2& pt,
                                     const Rat_point_2& q) const
  {
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
   * \return The shifted segemnt, represented as an x-monotone conic curve.
   */
  X_monotone_curve_2 _shift_segment (const Rat_point_2& ps,
                                     const Rat_point_2& pt,
                                     const Alg_point_2& q) const
  {
    Alg_point_2  shift_ps = Alg_point_2 (nt_traits.convert (ps.x()) + q.x(), 
                                         nt_traits.convert (ps.y()) + q.y());
    Alg_point_2  shift_pt = Alg_point_2 (nt_traits.convert (pt.x()) + q.x(), 
                                         nt_traits.convert (pt.y()) + q.y());

    return (X_monotone_curve_2 (shift_ps, shift_pt));
  }

  /*!
   */
  X_monotone_curve_2 _shift_arc (const Rat_circle_2& pc,
                                 const Alg_point_2& ps,
                                 const Alg_point_2& pt,
                                 const Rat_point_2& q) const
  {
  }

};


CGAL_END_NAMESPACE

#endif
