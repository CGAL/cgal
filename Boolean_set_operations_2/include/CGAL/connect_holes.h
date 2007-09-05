// Copyright (c) 2007  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
//
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>

#ifndef CGAL_CONNECT_HOLES_H
#define CGAL_CONNECT_HOLES_H

#include <CGAL/basic.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <list>
#include <set>

CGAL_BEGIN_NAMESPACE

template <class HANDLE>
struct _Less_handle
{
  bool operator() (const HANDLE& vh1, const HANDLE& vh2) const
  {
    return (&(*vh1) < &(*vh2));
  }
};

/*!
 * Connect the given polygon with holes, turning it into a sequence of
 * points, where the holes are connceted to the outer boundary using
 * zero-width passages.
 * For example:
 *              Input                             Output
 *  +----------------------------+    +-----*---------------*------+
 *  |                            |    |     |               |      |
 *  |     +------+        +--+   |    |     +------+        +--+   |
 *  |     |      |        |  |   |    |     |      |        |  |   |
 *  |     +------+         \ |   |    |     +----*-+         \ |   |
 *  |                       \|   |    |          |            \|   |
 *  |          +----+            |    |          +----+            |
 *  |         /     |            |    |         /     |            |
 *  |        +------+            |    |        +------+            |
 *  |                            |    |                            |
 *  +----------------------------+    +----------------------------+
 *
 * \param pwh The polygon with holes.
 * \param oi Output: An output iterator for the points.
 * \pre The polygons has an outer boundary.
 * \return A past-the-end iterator of the points.
 */
template <class Kernel, class Container, class OutputIterator>
OutputIterator connect_holes(const Polygon_with_holes_2<Kernel,
                             Container>& pwh,
                             OutputIterator oi)
{
  typedef Polygon_2<Kernel,Container>              Polygon_2;
  typedef Polygon_with_holes_2<Kernel,Container>   Polygon_with_holes_2;
  typedef Arr_segment_traits_2<Kernel>             Traits_2;
  typedef typename Kernel::Point_2                 Point_2;
  typedef typename Traits_2::X_monotone_curve_2    Segment_2;
  typedef Arrangement_2<Traits_2>                  Arrangement_2;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;

  CGAL_precondition (! pwh.is_unbounded());

  // In case the polygon has not holes, just go over its outer boundary
  // and report the points along it.
  const Polygon_2&                       outer_pgn = pwh.outer_boundary();
  typename Polygon_2::Vertex_circulator  first_v, curr_v, next_v;

  if (! pwh.has_holes())
  {
    first_v = curr_v = outer_pgn.vertices_circulator();
    
    do
    {
      *oi = *curr_v;
      ++oi;
      ++curr_v;

    } while (curr_v != first_v);

    return (oi);
  }

  // Go over the outer boundary of the polygon and store its edges.
  std::list<Segment_2>                   segments;

  first_v = curr_v = outer_pgn.vertices_circulator();
  do
  {
    next_v = curr_v;
    ++next_v;

    segments.push_back (Segment_2 (*curr_v, *next_v));
    curr_v = next_v;
    
  } while (curr_v != first_v);

  // Go over the holes and accumulate their edges as well.
  typename Polygon_with_holes_2::Hole_const_iterator   hole_it;

  for (hole_it = pwh.holes_begin(); hole_it != pwh.holes_end(); ++hole_it)
  {
    first_v = curr_v = hole_it->vertices_circulator();
    do
    {
      next_v = curr_v;
      ++next_v;

      segments.push_back (Segment_2 (*curr_v, *next_v));
      curr_v = next_v;

    } while (curr_v != first_v);
  }

  // Construct the arrangement of all segments.
  Arrangement_2         arr;

  insert (arr, segments.begin(), segments.end());

  // The resulting arrangment contains a single holes in the unbounded face,
  // which comprises a face f, with several holes in its interior.
  // Go over these holes and pick the topmost vertex in each hole.
  const Face_handle                      uf = arr.unbounded_face();
  typename Arrangement_2::Hole_iterator  f_hole_it = uf->holes_begin();
  const Face_handle                      f = (*f_hole_it)->twin()->face();
  typename Arrangement_2::Ccb_halfedge_circulator
                                         first, circ;
  Kernel                                 ker;
  typename Kernel::Compare_y_2           comp_y = ker.compare_y_2_object();
  typename Kernel::Compare_x_2           comp_x = ker.compare_x_2_object();
  Vertex_handle                          v_top;
  Comparison_result                      res;
  std::set<Vertex_const_handle,
           _Less_handle<Vertex_const_handle> >  top_vertices;

  for (f_hole_it = f->holes_begin(); f_hole_it != f->holes_end(); ++f_hole_it)
  {
    // Locate the topmost vertex in the current hole. In case of two (or more)
    // vertices with maximal y-coordinate, select the leftmost one.
    first = circ = *f_hole_it;
    v_top = circ->target();
    do
    {
      ++circ;

      res = comp_y (circ->target()->point(), v_top->point());

      if (res == CGAL::LARGER ||
          (res == CGAL::EQUAL &&
           comp_x (circ->target()->point(), v_top->point()) == CGAL::SMALLER))
      {
          v_top = circ->target();
      }
    } while (circ != first);

    top_vertices.insert (Vertex_const_handle (v_top));
  }

  // Perform "vertical ray shooting" from each arrangement vertex, locating
  // the features that lie below and above it.
  typedef std::list<std::pair<Vertex_const_handle,
                              std::pair<Object, Object> > >   Ray_shoot_list;

  Ray_shoot_list                     vrs_list;
  typename Ray_shoot_list::iterator  vrs_iter;

  decompose (arr, std::back_inserter (vrs_list));

  // Go over the results of the batched vertical ray-shooting query.
  Vertex_const_handle                v;
  Vertex_handle                      v_above;
  Halfedge_const_handle              he;
  Halfedge_handle                    he_above;
  typename Kernel::Direction_2          dir_up (0, 1);
  typename Kernel::Construct_ray_2      ray = ker.construct_ray_2_object();
  typename Kernel::Construct_segment_2  segment =
                                            ker.construct_segment_2_object();
  typename Kernel::Intersect_2          intersect = ker.intersect_2_object();
  Object                                obj;
  Point_2                               ip;
  bool                                  assign_success;

  for (vrs_iter = vrs_list.begin(); vrs_iter != vrs_list.end(); ++vrs_iter)
  {
    if (top_vertices.find (vrs_iter->first) == top_vertices.end())
      continue;

    v_top = arr.non_const_handle (vrs_iter->first);

    // In case the current vertex is a top vertex in one of the holes,
    // add a vertical segment connecting it with the feature above it.
    if (CGAL::assign (v, vrs_iter->second.second))
    {
      // v_top lies below a vertex v_above: Connect these two vertices.
      v_above = arr.non_const_handle (v);

      arr.insert_at_vertices (Segment_2 (v_top->point(), v_above->point()),
                              v_top, v_above);
    }
    else if (CGAL::assign (he, vrs_iter->second.second))
    {
      // v_top lies below the interior of the hafledge he_above:
      // Find the intersection of this halfegde with a vertical ray
      // emanating from v_top.
      he_above = arr.non_const_handle (he);

      obj = intersect (ray (v_top->point(), dir_up),
                       segment (he_above->source()->point(),
                                he_above->target()->point()));

      assign_success = CGAL::assign (ip, obj);
      CGAL_assertion (assign_success);

      if (assign_success)
      {
        // Split he_above at the computed intersection point.
        arr.split_edge (he_above,
                        Segment_2 (he_above->source()->point(), ip),
                        Segment_2 (ip, he_above->target()->point()));

        // Now he_above is split such that it becomes the predecessor
        // halfedge for the insertion of the vertical segment connecting
        // v_top and ip.
        arr.insert_at_vertices (Segment_2 (v_top->point(), ip),
                                he_above, v_top);
      }
    }
    else
    {
      // We should never reach here.
      CGAL_assertion_msg (false,
                          "top vertex is located in an unbounded face.");
    }
  }

  // The holes of the face f are now all connected to it outer boundary.
  // Go over this boundary and report the vertices along it.
  // Note that we start with a vertex located on the original outer boundary.
  first = f->outer_ccb();

  while (first->twin()->face() != uf)
    ++first;

  circ = first;
  do
  {
    *oi = circ->target()->point();
    ++oi;
    ++circ;

  } while (circ != first);

  return (oi);
}

CGAL_END_NAMESPACE

#endif
