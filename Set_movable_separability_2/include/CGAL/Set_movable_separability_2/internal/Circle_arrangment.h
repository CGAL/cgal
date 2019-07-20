// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s): Shahar    <shasha94@gmail.com>
//            Efi Fogel <efif@gmail.com>

#ifndef CGAL_SET_MOVABLE_SEPARABILITY_2_INTERNAL_CIRCLE_ARRANGMENT_H
#define CGAL_SET_MOVABLE_SEPARABILITY_2_INTERNAL_CIRCLE_ARRANGMENT_H

#include <CGAL/license/Set_movable_separability_2.h>


#include <CGAL/enum.h>

#include <iostream>
#include <list>

/* Legend:
 * point = Represented as Direction_2. It is the intersection between the
 *   fitting Direction_2 and the unit circle
 *
 * Arc = Represented as A pair of point. clockwise arc between the first
 *   point and the second point. (each of its sides might be open or closed)
 *
 * SegmentOuterCircle  = Arc that represent all the directions that points
 *   out from the polygon if it start from the
 *   fitting segment. This arc is always open half circle.
 */

/*! \Circle_arrangment
 * \brief This class represents an subdivision of the unit-circle into cells of
 * depth 0,1,2+ where depth is the number of inserted open half-circles inserted
 * that covers this cell in addition this class contains some static functions
 * that are in this class inorder of sharing its typedefs all the circle is
 * always covered by some cell. there can't be an hole
 */

namespace CGAL {
namespace Set_movable_separability_2 {
namespace internal {

template <typename Kernel>
class Circle_arrangment {
  typedef typename Kernel::Direction_2            Point;
  typedef std::pair<Point, Point>                 Arc;
  typedef typename CGAL::Polygon_2<Kernel>::Edge_const_iterator Edge_iter;

  /* Legend:
   * Point = Represented as Direction_2. It is the intersection between the
   *   fitting Direction_2 and the unit circle
   *
   * Arc = Represented as a pair of points. clockwise arc between the first
   *   point and the second point. (each of its sides might be open or closed)
   */

  /*! \fn bool is_open_direction_contained_in_arc(point p, bool is_counterclockwise, Arc A)
   * Checks whether an open epsilon area clockwise/counterclockwise from a point
   * p is contained in an arc s.
   * \param[in] p a point .
   * \param[in] is_counterclockwise true: we care about the counterclockwise
   *            epsilon area of p. false: same with clockwise
   * \param[in] A an Arc that should contain the epsilon area
   */
  bool is_open_direction_contained_in_arc(const Point p,
                                          const bool is_counterclockwise,
                                          const Arc A) const
  {
    if ((is_counterclockwise && (p == A.second)) ||
        (!is_counterclockwise && (p == A.first)))
      return false;
    auto cc_in_between = m_kernel.counterclockwise_in_between_2_object();
    return !cc_in_between(p, A.first, A.second);
  }

  /*! \fn bool is_a_contained_in_b(bool is_a_start_closed,bool is_a_end_closed, arc A, arc B)
 * \brief checks whether an arc A is contained in an arc B
 * \param[in] is_a_start_closed - do A contains its start point (clockwise)
 * \param[in] is_a_end_closed - do A contains its end point (clockwise)
 * \param[in] A - an arc
 * \param[in] B - an *open* arc
 */
  bool is_a_contained_in_b(const bool is_a_start_closed,
                           const bool is_a_end_closed,
                           const Arc A,const Arc B) const
  {
    //A is closed, B is open and they share an vertex -> A not contained in B
    if ((is_a_start_closed &&(A.first == B.first)) ||
        (is_a_end_closed && (A.second == B.second)))
      return false;
    if ((A.first == B.second) || (B.first == A.second)) return false;
    auto cc_in_between = m_kernel.counterclockwise_in_between_2_object();
    return (!cc_in_between(A.first, B.first, B.second) &&
            !cc_in_between(A.second, B.first, B.second) &&
            !cc_in_between(A.first, B.first, A.second));
  }

  /*! \Circle_arrangment_edge
   * This class represents a cells (a point or an arc) of depth 0,1,2+ in the
   * Circle_arrangment where depth the number of inserted open half-circles
   * inserted that cover this cell
   * This edge (cell) is described by the first point of the edge (clockwise).
   * The last point can be deduced by the next instance of
   * Circle_arrangment_edge in the list in Circle_arrangment
   * this class also keeps the cell depth.
   */
  class Circle_arrangment_edge {
  public:
    bool m_start_is_closed;

    Point m_edge_start_angle; // the end is the start of the next edge

    uint8_t m_count;          // no. of outer circles that cover the edge (0/1/2+)

    Edge_iter m_edge_iter;      // the iterator of the polygon edge the open
    // half-circle of which covers this cell.
    // only relevant if m_count ==1

    /*! \ctor Circle_arrangment_edge(point edge_start_angle, const Edge_iter edge_iter, bool start_is_closed,bool set_count_to_one=true)
     * Creates a new edge (Arc), this edge count must be 0 or 1
     * \param[in] edge_start_angle the first point of the arc (clockwise)
     * \param[in] edge_iter the iterator of the polygon edge who's open
     *            half-circle covers this cell - only relevant if m_count == 1
     * \param[in] start_is_closed - is the point edge_start_angle contained in
     *            this cell
     * \param[in] set_count_to_one to set the m_count to one (or zero if this
     *            var is false)
     */
    Circle_arrangment_edge(const Point edge_start_angle,
                           const Edge_iter edge_iter,
                           const bool start_is_closed,
                           const bool set_count_to_one = true)
    {
      this->m_start_is_closed = start_is_closed;
      this->m_edge_start_angle = edge_start_angle;
      this->m_count = (int) set_count_to_one;
      this->m_edge_iter = edge_iter;
    }

    /*! \fn void plusplus(const Edge_iter edge_ite)
     * Adds new polygon edge who's open half-circle covers this cell
     * \param[in] edge_ite - the iterator of this edge
     * increase the edge m_count by one (if it is 2+, it will stay 2+)
     * set this new edge to be the one covers the cell if the m_count was zero
     * before. (only relevant if now m_count == 1)
     */
    void plusplus(const Edge_iter edge_iter)
    {
      if (this->m_count ==0) {
        this->m_edge_iter = edge_iter;
        this->m_count = 1;
      }
      else if(this->m_count ==1) this->m_count = 2;
    }
    bool is_covered() { return m_count == 2; }
  };

  typedef typename std::list<Circle_arrangment_edge> Circle_edges;

  //! The kernel to use.
  const Kernel& m_kernel;

  Circle_edges m_edges;

  /*! \fn void insert_if_legal(const Circle_edge_iterator cur_it,
   *                           const Circle_edge_iterator next_it,
   *                           const Circle_arrangment_edge &edge)
   * Adds new edge to the arrangement if it won't create some empty edges
   * \param[in] cur_it iterator to the edge before where the new edge should be
   *        inserted
   * \param[in] next_it iterator to the edge after where the new edge should be
   *        inserted
   * \param[in] edge the new edge that should be inserted
   *
   * Notice that next_it is redundant since it can be deduced from cur_it.
   * But it was easier for me to just send it as well.
   */
  template <typename InputIterator>
  void insert_if_legal(const InputIterator cur_it,
                       InputIterator next_it,
                       const Circle_arrangment_edge& edge)
  {
    if (((edge.m_start_is_closed && !next_it->m_start_is_closed) ||
         (edge.m_edge_start_angle != next_it->m_edge_start_angle)) &&
        ((cur_it->m_start_is_closed && !edge.m_start_is_closed) ||
         (edge.m_edge_start_angle != cur_it->m_edge_start_angle)))
    {
      m_edges.insert(next_it, edge);
      return;
    }
  }

  /*! \fn void merge_adjacent_2_edges_and_remove_empty()
   * \brief merge all the arcs that are adjacent and of depth 2+
   * it doesn't merge the first and last ones since it is easier this way.
   */
  void merge_adjacent_2_edges_and_remove_empty()
  {
    bool in_two_edge(false);
    for (auto it = m_edges.begin(); it != m_edges.end();) {
      if (it->is_covered()) {
        if (in_two_edge) {
          it = m_edges.erase(it);
          continue;
        }
        in_two_edge = true;
      }
      else {
        in_two_edge = false;
      }
      ++it;
    }
  }

public:
  /*! \ctor Circle_arrangment(arc first_segment_outer_circle)
   * Creates an arrangement on circle with two edges the one covered by
   * first_segment_outer_circle and the other one
   * \param[in] first_segment_outer_circle the outer circle of the first segment
   *            of the polygon.
   * Notice that you might consider implementing the ctor as an full circle of
   * depth 0, but it was much easier for me to ignore the case where the all
   * circle is a single arc, so I choose this implementation.
   */
  Circle_arrangment(const Kernel& kernel, const Arc first_segment_outer_circle,
                    const Edge_iter edge_iter) :
    m_kernel(kernel)
  {
    m_edges.push_back(Circle_arrangment_edge(first_segment_outer_circle.first,
                                             edge_iter, false));
    m_edges.push_back(Circle_arrangment_edge(first_segment_outer_circle.second,
                                             edge_iter, true, false));
  }

  /*! \fn add_segment_outer_circle(arc segment_outer_circle, const Edge_iter edge_ite)
   * Updates the arrangement in respect to a new segment outer open circle
   * \param[in] segment_outer_circle - the outer circle of the current segment of
   *        the polygon.
   * \param[in] edge_ite this segment iterator
   * This is the main funtion of this code. It separates the cells in which the
   * endpoints of the new arc is contained to two parts and increase m_count
   * for all the cells that the new arc covers. In the end the function
   * merge_adjacent_2_edges_and_remove_empty is called to remove redundant cells
   */
  void add_segment_outer_circle(const Arc segment_outer_circle,
                                const Edge_iter edge_iter)
  {
    Arc edge;
    bool is_start_closed_segment = m_edges.begin()->m_start_is_closed;
    bool is_end_closed_segment ;
    edge.first = m_edges.begin()->m_edge_start_angle;
    //edge.second  ;
    auto next_it = m_edges.begin();
    auto it = m_edges.begin();
    bool done = false;
    while (!done) {
      it = next_it;
      next_it = it;
      ++next_it;
      if (next_it == m_edges.end()) {
        done = true;
        next_it = m_edges.begin();
      }

      is_start_closed_segment =it->m_start_is_closed;
      is_end_closed_segment = !next_it->m_start_is_closed;
      edge.first = it->m_edge_start_angle;
      edge.second =next_it->m_edge_start_angle;

      if (it->m_count == 2) continue;
      if (is_a_contained_in_b(is_start_closed_segment, is_end_closed_segment,
                              edge, segment_outer_circle))
      {
        it->plusplus(edge_iter);
        continue;
      }
      bool is_start_contained =
        is_open_direction_contained_in_arc(segment_outer_circle.first, true,
                                           edge);
      bool is_end_contained =
        is_open_direction_contained_in_arc(segment_outer_circle.second, false,
                                           edge);
      // o~~~~~~~~~~~~o  = new arc
      // ?------------?  = "old" arc (the edge from the array)
      if (is_start_contained) {
        if (is_end_contained) {
          auto cc_in_between = m_kernel.counterclockwise_in_between_2_object();
          bool isordered = !cc_in_between(segment_outer_circle.second,
                                          segment_outer_circle.first,
                                          edge.second);
          if (isordered) {
            //      o~~~~~~~~~~~~o
            // ?-----------------------?
            // __________________________
            // ?----c
            //      o~-~-~-~-~-~-o
            //                   c-----?
            Circle_arrangment_edge edge2 = *it;
            edge2.m_start_is_closed = false;
            edge2.m_edge_start_angle = segment_outer_circle.first;
            edge2.plusplus(edge_iter);
            this->insert_if_legal(it, next_it, edge2);
            Circle_arrangment_edge edge3 = *it;
            edge3.m_start_is_closed = true;
            edge3.m_edge_start_angle = segment_outer_circle.second;
            this->insert_if_legal(it,next_it,edge3);
          }
          else {
            // ...~~~~~~~~~o   o~~~~~~~~~~... (round)
            //         ?-----------?
            // __________________________
            //         ?~-~o
            //             c---c
            //                 o-~-?
            Circle_arrangment_edge edge2 = *it;
            edge2.m_start_is_closed = true;
            edge2.m_edge_start_angle = segment_outer_circle.second;
            this->insert_if_legal(it, next_it, edge2);
            Circle_arrangment_edge edge3 = *it;
            edge3.m_start_is_closed = false;
            edge3.m_edge_start_angle = segment_outer_circle.first;
            edge3.plusplus(edge_iter);
            this->insert_if_legal(it, next_it, edge3);
            it->plusplus(edge_iter);
          }
        }
        else {
          //      o~~~~~~~~~~~~o
          // ?-----------?
          //_____________________
          // ?----c
          //      o-~-~-~?
          Circle_arrangment_edge edge2 = *it;
          edge2.m_start_is_closed = false;
          edge2.m_edge_start_angle = segment_outer_circle.first;
          edge2.plusplus(edge_iter);
          this->insert_if_legal(it, next_it, edge2);
        }
      }
      else {
        if (is_end_contained) {
          // o~~~~~~~~~~~~o
          //        ?------------?
          //_____________________
          //        ?-~-~-~-o
          //                c----?
          Circle_arrangment_edge edge2 = *it;
          edge2.m_start_is_closed = true;
          edge2.m_edge_start_angle = segment_outer_circle.second;
          it->plusplus(edge_iter);
          this->insert_if_legal(it, next_it, edge2);
        }
        //else -  no intersection, do noting
      }
    }
    merge_adjacent_2_edges_and_remove_empty();
  }

#if 0
  // debug function
  void printArrangement()
  {
    for (auto it = m_edges.begin(); it != m_edges.end(); ++it) {
      if (it->m_start_is_closed) std::cout<<")[";
      else std::cout << "](";
      std::cout << it->m_edge_start_angle;
      std::cout << ","<<(int)it->m_count;
    }
    std::cout << "\n\n";
  }
#endif

  /*! \fn void get_all_1_edges(OutputIterator oi)
   * Insert to oi all the cells in depth 1 i.e. the cells that represent legal
   * pullout directions
   * \param[in, out] oi the output iterator to put the cells in
   * Puts in oi var of type pair<size_t, std::pair<Kernel::Direction_2,
   * Kernel::Direction_2 > > foreach valid top edge.
   * Should only be called after all of the polygon edges where inserted.
   */
  template <typename OutputIterator>
  OutputIterator get_all_1_edges(OutputIterator oi)
  {
    for (auto it = m_edges.begin(); it != m_edges.end();) {
      if ((*it).m_count == 1) {
        std::pair<Edge_iter, Arc> edge;
        edge.first = (*it).m_edge_iter;
        edge.second.first = (*it).m_edge_start_angle;
        ++it;
        edge.second.second =
          (*((it == m_edges.end()) ? m_edges.begin() : (it))).m_edge_start_angle;
        *oi++ = edge;
      }
      else
      {
        ++it;
      }
    }
    return oi;
  }

  /*! \fn bool all_is_covered_twice()
   * Before running this run merge_adjacent_2_edges_and_remove_empty() or
   * add_segment_outer_circle() which calls
   * merge_adjacent_2_edges_and_remove_empty().
   *
   * The funtions checks that the whole circle is a single cell, which can
   * happen only if this cell is of depth 2, so there is no need to check the
   * depth as well.
   * \return if all of the arrangement is in depth 2+
   */
  bool all_is_covered_twice() { return m_edges.size() == 1; }
};

} // namespace internal
} // namespace Set_movable_separability_2
} // namespace CGAL

#endif
