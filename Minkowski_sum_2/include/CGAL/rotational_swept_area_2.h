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

#ifndef CGAL_ROTATIONAL_SWEPT_AREA_2_H
#define CGAL_ROTATIONAL_SWEPT_AREA_2_H

#include <CGAL/basic.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <list>

CGAL_BEGIN_NAMESPACE

/*!
 * Compute the area swept by a simple linear polygon when rotated in a
 * counterclockwise direction from orientation theta1 to orientation theta2
 * about a given point.
 * \param pgn The polygon.
 * \param p The center of rotation.
 * \param sin_theta1, cos_theta1 The sine and cosine of the source orientation.
 * \param sin_theta2, cos_theta2 The sine and cosine of the target orientation.
 */
template <class Kernel, class Container>
typename Gps_circle_segment_traits_2<Kernel>::Polygon_2
rotational_swept_area_2 (const Polygon_2<Kernel,Container>& pgn,
                         const typename Kernel::Point_2& p,
                         const typename Kernel::FT& sin_theta1,
                         const typename Kernel::FT& cos_theta1,
                         const typename Kernel::FT& sin_theta2,
                         const typename Kernel::FT& cos_theta2)
{
  typedef typename Kernel::FT                         NT;
  typedef typename Kernel::Point_2                    Point_2;
  typedef Polygon_2<Kernel,Container>                 Polygon_2;
  typedef typename Polygon_2::Vertex_const_circulator Vertex_circulator;

  typedef Gps_circle_segment_traits_2<Kernel>         Traits_2;
  typedef typename Traits_2::Point_2                  Tr_point_2;
  typedef typename Traits_2::Curve_2                  Curve_2;
  typedef typename Traits_2::X_monotone_curve_2       X_monotone_curve_2;
  typedef typename Traits_2::Polygon_2                Tr_polygon_2;

  typedef Arrangement_2<Traits_2>                     Arrangement_2;

  // Make sure that we deal with valid rotations.
  CGAL_precondition (CGAL::square(sin_theta1) + 
                     CGAL::square(cos_theta1) == NT(1));

  CGAL_precondition (CGAL::square(sin_theta2) + 
                     CGAL::square(cos_theta2) == NT(1));

  // Prepare circulators over the polygon vertices.
  const bool            forward = (pgn.orientation() == COUNTERCLOCKWISE);
  Vertex_circulator     first, curr, next;
  bool                  first_vertex = true;
  const NT              px = p.x(), py = p.y();
  Point_2               s1, t1;
  Point_2               s2, t2;
  std::list<Curve_2>    curves;

  first = pgn.vertices_circulator();
  curr = first; 
  next = first;
  
  // Go over the polygon vertices.
  do
  {
    // Get a circulator for the next vertex (in counterclockwise orientation).
    if (forward)
      ++next;
    else
      --next;
  
    // Let us denote the current edge by st. We compute s1 t1, the endpoints
    // rotated by theta1, and s2 t2, the endpoints rotated by theta2.
    if (first_vertex)
    {
      s1 = Point_2 (px + cos_theta1 * (curr->x() - px) - 
                         sin_theta1 * (curr->y() - py),
                    py + cos_theta1 * (curr->y() - py) +
                         sin_theta1 * (curr->x() - px));
      s2 = Point_2 (px + cos_theta2 * (curr->x() - px) - 
                         sin_theta2 * (curr->y() - py),
                    py + cos_theta2 * (curr->y() - py) + 
                         sin_theta2 * (curr->x() - px));

      first_vertex = false;
    }

    t1 = Point_2 (px + cos_theta1 * (next->x() - px) - 
                       sin_theta1 * (next->y() - py),
                  py + cos_theta1 * (next->y() - py) +
                       sin_theta1 * (next->x() - px));
    t2 = Point_2 (px + cos_theta2 * (next->x() - px) - 
                       sin_theta2 * (next->y() - py),
                  py + cos_theta2 * (next->y() - py) + 
                       sin_theta2 * (next->x() - px));

    // Construct the segments s1 t1 and s2 t2.
    curves.push_back (Curve_2 (s1, t1));
    curves.push_back (Curve_2 (s2, t2));

    // Construct the ciruclar arc centered at p and connecting s1 and s2.
    const NT                    r_sqr = CGAL::square (curr->x() - px) + 
                                        CGAL::square (curr->y() - py);
    typename Kernel::Circle_2   circ (p, r_sqr, COUNTERCLOCKWISE);

    Tr_point_2   ps1 (s1.x(), s1.y());
    Tr_point_2   ps2 (s2.x(), s2.y());

    curves.push_back (Curve_2 (circ, ps1, ps2));

    // Proceed to the next vertex.
    curr = next;
    s1 = t1;
    s2 = t2;
  
  } while (curr != first);

  // Construct the arrangement of all curves.
  Arrangement_2         arr;

  insert_curves (arr, curves.begin(), curves.end());

  // The resulting arrangement should contain a single hole in the unbounded
  // face. We return the x-monotone curves that form the boundary of this
  // hole. 
  typename Arrangement_2::Face_const_handle       uf = arr.unbounded_face();
  typename Arrangement_2::Hole_const_iterator     hole_it = uf->holes_begin();
  typename Arrangement_2::Ccb_halfedge_const_circulator  first_circ, circ;
  std::list<X_monotone_curve_2>                          boundary_curves;

  circ = first_circ = *hole_it;
  do
  {
    boundary_curves.push_back (circ->curve());
    ++circ;

  } while (circ != first_circ);
    
  // Make sure that there is a single hole in the unbounded face.
  ++hole_it;
  CGAL_assertion (hole_it == uf->holes_end());

  return (Tr_polygon_2 (boundary_curves.rbegin(), 
                        boundary_curves.rend()));
}

CGAL_END_NAMESPACE

#endif
