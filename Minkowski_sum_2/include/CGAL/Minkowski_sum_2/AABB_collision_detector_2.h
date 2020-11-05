// Copyright (c) 2015  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Sebastian Morr    <sebastian@morr.cc>

#ifndef CGAL_AABB_COLLISION_DETECTOR_2_H
#define CGAL_AABB_COLLISION_DETECTOR_2_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Minkowski_sum_2/AABB_tree_with_join.h>
#include <CGAL/Minkowski_sum_2/AABB_traits_2.h>
#include <CGAL/Minkowski_sum_2/AABB_segment_2_primitive.h>

namespace CGAL {

// Tests whether two polygons P and Q overlap for different translations of Q.
template <class Kernel_, class Container_>
class AABB_collision_detector_2
{
public:
  typedef Kernel_                                       Kernel;
  typedef Container_                                    Container;

  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Vector_2                     Vector_2;
  typedef typename CGAL::Polygon_2<Kernel, Container>   Polygon_2;
  typedef typename CGAL::Polygon_with_holes_2<Kernel, Container>
                                                        Polygon_with_holes_2;
  typedef typename Polygon_2::Edge_const_iterator       Edge_iterator;
  typedef AABB_segment_2_primitive<Kernel, Edge_iterator, Polygon_with_holes_2>
                                                        Tree_segment_2;
  typedef AABB_traits_2<Kernel, Tree_segment_2>         Tree_traits;
  typedef AABB_tree_with_join<Tree_traits>              Tree_2;

public:
  AABB_collision_detector_2(const Polygon_with_holes_2& p,
                            const Polygon_with_holes_2& q) :
    m_p(q), m_q(p)
  {
    m_stationary_tree.insert(p.outer_boundary().edges_begin(),
                             p.outer_boundary().edges_end());

    typename Polygon_with_holes_2::Hole_const_iterator it = p.holes_begin();
    while (it != p.holes_end())
    {
      m_stationary_tree.insert(it->edges_begin(), it->edges_end());
      ++it;
    }

    m_translating_tree.insert(q.outer_boundary().edges_begin(),
                              q.outer_boundary().edges_end());

    it = q.holes_begin();
    while (it != q.holes_end())
    {
      m_translating_tree.insert(it->edges_begin(), it->edges_end());
      ++it;
    }
  }

  // Returns true iff the polygons' boundaries intersect or one polygon is
  // completely inside of the other one. Q is translated by t.
  bool check_collision(const Point_2 &t)
  {
    if (m_stationary_tree.do_intersect(m_translating_tree, t)) return true;

    // If t_q is inside of P, or t_p is inside of Q, one polygon is completely
    // inside of the other.

    // Obtain a point on the boundary of m_q:
    Point_2 t_q = (! m_q.outer_boundary().is_empty()) ?
      *m_q.outer_boundary().vertices_begin() - Vector_2(ORIGIN, t) :
      *m_q.holes_begin()->vertices_begin() - Vector_2(ORIGIN, t);

    // Obtain a point on the boundary of m_p:
    Point_2 t_p = (! m_p.outer_boundary().is_empty()) ?
      *m_p.outer_boundary().vertices_begin() + Vector_2(ORIGIN, t) :
      *m_p.holes_begin()->vertices_begin() + Vector_2(ORIGIN, t);

    // Use bounded_side_2() instead of on_bounded_side() because the latter
    // checks vor simplicity every time.
    bool in_mp(true);
    if (! m_p.outer_boundary().is_empty())
      in_mp =
        bounded_side_2(m_p.outer_boundary().vertices_begin(),
                       m_p.outer_boundary().vertices_end(), t_q,
                       m_p.outer_boundary().traits_member()) == ON_BOUNDED_SIDE;
    if (m_p.number_of_holes() == 0) {
      if (in_mp) return true;
    }
    bool in_mq(true);
    if (! m_q.outer_boundary().is_empty())
      in_mq =
        bounded_side_2(m_q.outer_boundary().vertices_begin(),
                       m_q.outer_boundary().vertices_end(), t_p,
                       m_q.outer_boundary().traits_member()) == ON_BOUNDED_SIDE;
    if (m_q.number_of_holes() == 0) {
      if (in_mq) return true;
    }
    if (!in_mq && !in_mp) return false;
    if (in_mp) {
      for (typename Polygon_with_holes_2::Hole_const_iterator it =
             m_p.holes_begin(); it != m_p.holes_end(); ++it)
      {
        if (bounded_side_2(it->vertices_begin(), it->vertices_end(), t_q,
                           it->traits_member()) == ON_BOUNDED_SIDE)
          return false;
      }
      return true;
    }
    for (typename Polygon_with_holes_2::Hole_const_iterator it =
           m_q.holes_begin(); it != m_q.holes_end(); ++it)
    {
      if (bounded_side_2(it->vertices_begin(), it->vertices_end(), t_p,
                         it->traits_member()) == ON_BOUNDED_SIDE)
        return false;
    }
    return true;
  }

private:
  Tree_2 m_stationary_tree;
  Tree_2 m_translating_tree;
  const Polygon_with_holes_2& m_p;
  const Polygon_with_holes_2& m_q;
};

} // namespace CGAL

#endif
