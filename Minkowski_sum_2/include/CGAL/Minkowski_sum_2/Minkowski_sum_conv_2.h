// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// Author(s) : Ron Wein   <wein_r@yahoo.com>
//             Efi Fogel  <efifogel@gmail.com>

#ifndef CGAL_MINKOWSKI_SUM_CONV_H
#define CGAL_MINKOWSKI_SUM_CONV_H

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Minkowski_sum_2/Labels.h>
#include <CGAL/Minkowski_sum_2/Arr_labeled_traits_2.h>
#include <CGAL/Minkowski_sum_2/Union_of_segment_cycles_2.h>

namespace CGAL {

/*! \class
 * A class for computing the Minkowski sum of two simple polygons based on the
 * convolution of their boundaries.
 */
template <typename Kernel_, typename Container_>
class Minkowski_sum_by_convolution_2 {
public:
  typedef Kernel_                                        Kernel;
  typedef CGAL::Polygon_2<Kernel, Container_>            Polygon_2;

private:
  // Kernel types:
  typedef typename Kernel::Point_2                       Point_2;
  typedef typename Kernel::Vector_2                      Vector_2;
  typedef typename Kernel::Direction_2                   Direction_2;

  // Kernel functors:
  typedef typename Kernel::Equal_2                       Equal_2;
  typedef typename Kernel::Construct_translated_point_2  Translate_point_2;
  typedef typename Kernel::Construct_vector_2            Construct_vector_2;
  typedef typename Kernel::Construct_direction_2         Construct_direction_2;
  typedef typename Kernel::Construct_opposite_line_2     Opposite_line_2;
  typedef typename Kernel::Orientation_2                 Compute_orientation_2;
  typedef typename Kernel::Compare_xy_2                  Compare_xy_2;
  typedef typename Kernel::Counterclockwise_in_between_2 Ccw_in_between_2;

  // Polygon-related types:
  typedef typename Polygon_2::Vertex_circulator          Vertex_circulator;
  typedef std::pair<Vertex_circulator, unsigned int>     Vertex_ref;
  typedef std::pair<Vertex_ref, Vertex_ref>              Anchor;
  typedef std::list<Anchor>                              Anchors_queue;

  /*! \class
   * An auxiliary class for representing labels of convolved vertex pairs.
   */
  class Convolution_label {
  private:
    unsigned int index1;       // Vertex index of the first polygon.
    unsigned int index2;       // Vertex index of the second polygon.
    unsigned int move_on;      // On which polygon do we move.

  public:
    /*! Default constructor. */
    Convolution_label() : move_on(0) {}

    /*! Constructor with parameters. */
    Convolution_label(unsigned int ind1, unsigned int ind2,
                      unsigned int move) :
      index1(ind1),
      index2(ind2),
      move_on(move)
    {
      CGAL_precondition(move_on == 1 || move_on == 2);
    }

    /*! Less-then operator (for the usage of std::set). */
    bool operator<(const Convolution_label& label) const
    {
      if (index1 < label.index1) return (true);
      else if (index1 > label.index1) return (false);

      if (index2 < label.index2) return (true);
      else if (index2 > label.index2) return (false);

      return (move_on < label.move_on);
    }
  };

  typedef std::set<Convolution_label>                     Labels_set;

  // Traits-related types:
  typedef Arr_segment_traits_2<Kernel>                    Segment_traits_2;
  typedef Arr_labeled_traits_2<Segment_traits_2>          Traits_2;

  typedef typename Segment_traits_2::X_monotone_curve_2   Segment_2;
  typedef typename Traits_2::X_monotone_curve_2           Labeled_segment_2;
  typedef std::list<Labeled_segment_2>                    Segments_list;

  typedef Union_of_segment_cycles_2<Traits_2, Polygon_2>  Union_2;

  const Kernel* m_kernel;
  bool m_own_kernel;    // inidicates whether the kernel should be freed up.

  // Data members:
  Equal_2                 f_equal;
  Translate_point_2       f_add;
  Construct_vector_2      f_vector;
  Construct_direction_2   f_direction;
  Opposite_line_2         f_opp_line;
  Compute_orientation_2   f_orientation;
  Compare_xy_2            f_compare_xy;
  Ccw_in_between_2        f_ccw_in_between;

public:
  /*! Default constructor. */
  Minkowski_sum_by_convolution_2() :
    m_kernel(new Kernel),
    m_own_kernel(true)
  { init(); }

  /*! Constructor. */
  Minkowski_sum_by_convolution_2(const Kernel& kernel) :
    m_kernel(&kernel),
    m_own_kernel(false)
  { init(); }

  /*! Destructor */
  ~Minkowski_sum_by_convolution_2()
  {
    if (m_own_kernel) {
      if (m_kernel != NULL) {
        delete m_kernel;
        m_kernel = NULL;
      }
      m_own_kernel = false;
    }
  }

  /*! Initialize. */
  void init()
  {
    // Obtain kernel functors.
    f_equal = m_kernel->equal_2_object();
    f_add = m_kernel->construct_translated_point_2_object();
    f_vector = m_kernel->construct_vector_2_object();
    f_direction = m_kernel->construct_direction_2_object();
    f_opp_line = m_kernel->construct_opposite_line_2_object();
    f_orientation = m_kernel->orientation_2_object();
    f_compare_xy = m_kernel->compare_xy_2_object();
    f_ccw_in_between = m_kernel->counterclockwise_in_between_2_object();
  }

  /*!
   * Obtain the kernel
   * \return the kernel
   */
  const Kernel* kernel() const { return m_kernel; }

  /*!
   * Compute the Minkowski sum of two simple polygons.
   * Note that as the input polygons may not be convex, the Minkowski sum may
   * not be a simple polygon. The result is therefore represented as
   * the outer boundary of the Minkowski sum (which is always a simple polygon)
   * and a container of simple polygons, representing the holes inside this
   * polygon.
   * \param pgn1 The first polygon.
   * \param pgn2 The second polygon.
   * \param sum_bound Output: A polygon respresenting the outer boundary
   *                          of the Minkowski sum.
   * \param sum_holes Output: An output iterator for the holes in the sum,
   *                          represented as simple polygons.
   * \pre Both input polygons are simple.
   * \return A past-the-end iterator for the holes in the sum.
   */
  template <typename OutputIterator>
  OutputIterator operator()(const Polygon_2& pgn1, const Polygon_2& pgn2,
                            Polygon_2& sum_bound, OutputIterator sum_holes)
    const
  {
    CGAL_precondition(pgn1.is_simple());
    CGAL_precondition(pgn2.is_simple());

    // Prepare the vector of directions for the first polygon.
    const unsigned int n1 = static_cast<unsigned int>(pgn1.size());
    const bool forward1 = (pgn1.orientation() == COUNTERCLOCKWISE);
    std::vector<Direction_2> dirs1(n1);
    Vertex_circulator curr1, next1;
    unsigned int k1;

    next1 = curr1 = pgn1.vertices_circulator();
    for (k1 = 0; k1 < n1; ++k1) {
      if (forward1) ++next1;
      else --next1;

      dirs1[k1] = f_direction(f_vector(*curr1, *next1));
      curr1 = next1;
    }

    // Prepare the vector of directions for the second polygon.
    // Also prepare a list containing all reflex vertices of this polygon.
    const unsigned int n2 = static_cast<unsigned int>(pgn2.size());
    const bool forward2 = (pgn2.orientation() == COUNTERCLOCKWISE);
    std::vector<Direction_2> dirs2(n2);
    Vertex_circulator prev2, curr2, next2;
    Vertex_ref bottom_left;
    // bool                      is_convex2 = true;
    std::list<Vertex_ref> reflex_vertices;
    unsigned int k2;

    prev2 = next2 = curr2 = pgn2.vertices_circulator();

    if (forward2) --prev2;
    else ++prev2;

    for (k2 = 0; k2 < n2; ++k2) {
      if (forward2) ++next2;
      else --next2;

      if (k2 == 0)
        bottom_left = Vertex_ref(curr2, 0);
      else if (f_compare_xy(*curr2, *(bottom_left.first)) == SMALLER)
        bottom_left = Vertex_ref(curr2, k2);

      if (f_orientation(*prev2, *curr2, *next2) == RIGHT_TURN) {
        // We found a reflex vertex.
        // is_convex2 = false;
        reflex_vertices.push_back (Vertex_ref (curr2, k2));
      }

      dirs2[k2] = f_direction(f_vector(*curr2, *next2));
      prev2 = curr2;
      curr2 = next2;
    }

    // Add the bottom-left vertex of the second polygon to the reflex vertices.
    typename std::list<Vertex_ref>::iterator  reflex_it;

    reflex_vertices.push_front(bottom_left);

    // Construct the segments of the convolution cycles.
    unsigned int curr_id = 0;
    unsigned int cycles = 0;
    Segments_list conv_segments;
    Segments_list cycle;
    Labels_set used_labels;
    Anchor anchor;
    Anchors_queue queue;
    unsigned int loops;

    for (reflex_it = reflex_vertices.begin();
         reflex_it != reflex_vertices.end(); ++reflex_it)
    {
      // Get the current reflex vertex (or the bottom-left vertex).
      next2 = curr2 = reflex_it->first;
      k2 = reflex_it->second;

      if (forward2) ++next2;
      else --next2;

      // Search the first polygon for a vertex that contains (k1, k2, 1) in
      // a convolution cycle.
      next1 = curr1 = pgn1.vertices_circulator();
      for (k1 = 0; k1 < n1; ++k1) {
        if (forward1) ++next1;
        else --next1;

        if ((used_labels.count(Convolution_label(k1, k2, 1)) == 0 &&
             (f_ccw_in_between(dirs1[k1], dirs2[(n2 + k2 - 1) % n2],
                               dirs2[k2]) ||
              f_equal(dirs1[k1], dirs2[k2]))) ||
            (used_labels.count(Convolution_label(k1, k2, 2)) == 0 &&
             f_ccw_in_between(dirs2[k2], dirs1[(n1 + k1 - 1) % n1],
                              dirs1[k1])))
        {
          // Construct the current convolution cycle.
          queue.clear();
          queue.push_back(Anchor(Vertex_ref(curr1, k1),
                                 Vertex_ref(curr2, k2)));
          loops = 0;

          while (! queue.empty()) {
            // Pop the first pair of anchor vertices from the queue.
            anchor = queue.front();
            queue.pop_front();
            loops++;

            const Vertex_ref& vert1 = anchor.first;
            const Vertex_ref& vert2 = anchor.second;

            if (loops > 0 &&
                used_labels.count(Convolution_label(vert1.second,
                                                    vert2.second, 2)) != 0)
            {
              loops--;
              continue;
            }

            // Add a loop to the current convolution cycle.
            curr_id++;
            _convolution_cycle(curr_id,
                               n1, forward1, dirs1, vert1.first, vert1.second,
                               n2, forward2, dirs2, vert2.first, vert2.second,
                               used_labels, queue,
                               cycle);

            // Catenate the segments of the current loop to the convolution
            // list.
            if (cycle.empty()) {
              --loops;
            }
            else {
              conv_segments.splice(conv_segments.end(), cycle);
              CGAL_assertion(cycle.empty());
            }
          }
          ++cycles;
        }

        curr1 = next1;
      }
    }

    // Compute the union of the cycles that represent the Minkowski sum.
    Union_2     unite;

    sum_holes = unite(conv_segments.begin(), conv_segments.end(),
                      sum_bound, sum_holes);

    return (sum_holes);
  }

private:
  /*!
   * Compute a convolution cycle starting from tow given vertices.
   * \param cycle_id The index of the current cycle.
   * \param n1 The size of the first polygon.
   * \param forward1 Whether we move forward or backward on this polygon.
   * \param dirs1 The directions of the edges in the first polygon.
   * \param curr1 Points to the current vertex in the first polygon.
   * \param k1 The index of this vertex (between 0 and n1 - 1).
   * \param n2 The size of the second polygon.
   * \param forward2 Whether we move forward or backward on this polygon.
   * \param dirs2 The directions of the edges in the second polygon.
   * \param curr2 Points to the current vertex in the second polygon.
   * \param k2 The index of this vertex (between 0 and n2 - 1).
   * \param used_labels Input/output: The segment labels used so far.
   * \param queue A queue of anchor vertices for loops in the cycle.
   * \param cycle Output: An list of labeled segments that constitute the
   *                      convolution cycle.
   * \return A past-the-end iterator for the output segments.
   */
  void _convolution_cycle(unsigned int cycle_id,
                          unsigned int n1, bool forward1,
                          const std::vector<Direction_2>& dirs1,
                          Vertex_circulator curr1, unsigned int k1,
                          unsigned int n2, bool forward2,
                          const std::vector<Direction_2>& dirs2,
                          Vertex_circulator curr2, unsigned int k2,
                          Labels_set& used_labels,
                          Anchors_queue& queue,
                          Segments_list& cycle) const

  {
    // Update the circulator pointing to the next vertices in both polygons.
    unsigned int seg_index = 0;
    const unsigned int first1 = k1;
    const unsigned int first2 = k2;
    Vertex_circulator next1 = curr1;
    Vertex_circulator next2 = curr2;
    Convolution_label label;
    Point_2 first_pt, curr_pt, next_pt;
    bool inc1, inc2;
    Comparison_result res;
    const bool MOVE_ON_1 = true, MOVE_ON_2 = false;

    if (forward1) ++next1;
    else --next1;

    if (forward2) ++next2;
    else --next2;

    // Start constructing the convolution cycle from *curr1 and *curr2.
    curr_pt = first_pt = f_add(*curr1, f_vector(CGAL::ORIGIN, *curr2));
    do {
      // Determine on which polygons we should move.
      inc1 = false;
      inc2 = false;

      if (f_ccw_in_between(dirs1[k1], dirs2[(n2 + k2 - 1) % n2], dirs2[k2])) {
        inc1 = true;

        label = Convolution_label(k1, k2, 1);
        if (used_labels.count(label) != 0)
          inc1 = false;
      }

      if (f_ccw_in_between(dirs2[k2], dirs1[(n1 + k1 - 1) % n1], dirs1[k1])) {
        if (inc1) {
          // If we are about to increment the first polygon, add an anchor
          // to the queue. Next time we reach here we will increment the
          // second polygon (and proceed until reaching this point again and
          // closing the loop).
          label = Convolution_label(k1, k2, 2);
          if (used_labels.count(label) == 0) {
            queue.push_back(Anchor(Vertex_ref(curr1, k1),
                                   Vertex_ref(curr2, k2)));
          }
        }
        else {
          inc2 = true;

          label = Convolution_label(k1, k2, 2);
          if (used_labels.count(label) != 0)
            inc2 = false;
        }
      }

      if (! inc1 && ! inc2 && f_equal(dirs1[k1], dirs2[k2])) {
        inc1 = true;
        inc2 = true;

        label = Convolution_label(k1, k2, 1);
        if (used_labels.count(label) != 0) inc1 = false;

        if (inc1) label = Convolution_label((k1 + 1) % n1, k2, 2);
        else label = Convolution_label(k1, k2, 2);

        if (used_labels.count(label) != 0) inc2 = false;
      }

      CGAL_assertion(inc1 || inc2);

      // Act according to the increment flags.
      if (inc1) {
        // Translate the current edge of the first polygon to *curr2.
        next_pt = f_add(*next1, f_vector(CGAL::ORIGIN, *curr2));

        res = f_compare_xy(curr_pt, next_pt);
        CGAL_assertion(res != EQUAL);

        cycle.push_back(Labeled_segment_2(Segment_2(curr_pt, next_pt),
                                          X_curve_label((res == SMALLER),
                                                        cycle_id,
                                                        seg_index,
                                                        MOVE_ON_1)));
        used_labels.insert(Convolution_label(k1, k2, 1));
        ++seg_index;

        // Proceed to the next vertex of the first polygon.
        curr1 = next1;
        k1 = (k1 + 1) % n1;

        if (forward1) ++next1;
        else --next1;

        curr_pt = next_pt;
      }

      if (inc2) {
        // Translate the current edge of the second polygon to *curr1.
        next_pt = f_add(*next2, f_vector(CGAL::ORIGIN, *curr1));

        res = f_compare_xy(curr_pt, next_pt);
        CGAL_assertion(res != EQUAL);

        cycle.push_back(Labeled_segment_2(Segment_2(curr_pt, next_pt),
                                          X_curve_label((res == SMALLER),
                                                        cycle_id,
                                                        seg_index,
                                                        MOVE_ON_2)));
        used_labels.insert(Convolution_label(k1, k2, 2));
        seg_index++;

        // Proceed to the next vertex of the second polygon.
        curr2 = next2;
        k2 = (k2 + 1) % n2;

        if (forward2) ++next2;
        else --next2;

        curr_pt = next_pt;
      }

    } while ((k1 != first1) || (k2 != first2));

    CGAL_assertion(f_equal(curr_pt, first_pt));

    // Before moving un-necessary sub-cycles from the segment list, make sure
    // the list contains no "cyclic" sub-cylces. We do that by making sure that
    // the first and last segments of the list correspond to traversals of
    // different polygons.
    int steps = static_cast<int>(cycle.size());

    while ((steps > 0) &&
           cycle.front().label().get_flag() == cycle.back().label().get_flag())
    {
      cycle.push_back(cycle.front());
      cycle.pop_front();
      steps--;
    }

    if (steps == 0) {
      cycle.clear();
      return;
    }

    // Reduce un-necessary sub-cycles. This is done by scanning the segments
    // list for subsequences sharing a common move_on indicator. When we
    // encounter such a subsequence that equals the size of the corresponding
    // polygon, we can safely remove it from the convolution cycle.
    typename std::list<Labeled_segment_2>::iterator  first, curr;
    bool move_on;
    unsigned int count = 1;
    bool reduced_cycle = false;

    curr = first = cycle.begin();
    move_on = first->label().get_flag();

    curr->label().set_flag(false);
    ++curr;
    while (curr != cycle.end()) {
      if ((move_on == MOVE_ON_1 && count == n1) ||
          (move_on == MOVE_ON_2 && count == n2))
      {
        // We have discovered a sequence of moves on one of the polygon that
        // equals the polygon size, so we can remove this sequence.
        cycle.erase(first, curr);
        reduced_cycle = true;
        first = curr;
        move_on = first->label().get_flag();
        count = 1;
      }
      else {
        if (move_on == curr->label().get_flag()) {
          ++count;
        }
        else {
          first = curr;
          move_on = first->label().get_flag();
          count = 1;
        }
      }

      curr->label().set_flag(false);
      ++curr;
    }

    if ((move_on == MOVE_ON_1 && count == n1) ||
        (move_on == MOVE_ON_2 && count == n2))
    {
      cycle.erase(first, curr);
      reduced_cycle = true;
    }

    // Eliminate "antenna"s that occur in the cycle.
    typename std::list<Labeled_segment_2>::iterator next, after_next;
    bool is_last = false;

    next = curr = cycle.begin();
    ++next;

    do {
      // If the current segment is the last one in the cycle, its next segment
      // (in cyclic order) is the first one.
      if (next == cycle.end()) {
        is_last = true;
        next = cycle.begin();
        after_next = cycle.end();
      }
      else {
        after_next = next;
        ++after_next;
      }

      // Check whether the current segment and the next segment form an
      // "antenna". This can happen only if their supporting lines are
      // opposite (as we know curr's target equals next's source):
      if (f_equal(curr->line(), f_opp_line(next->line()))) {
        // In case curr's source equals next's target, then the two segments
        // are opposite and cancel one another. We can therefore remove them
        // both. Otherwise, we replace them with the segment directed from
        // curr's source to next's target.
        curr_pt = curr->source();
        next_pt = next->target();

        res = f_compare_xy(curr_pt, next_pt);

        if (res != EQUAL) {
          cycle.insert(curr,
                       Labeled_segment_2(Segment_2(curr_pt, next_pt),
                                         X_curve_label((res == SMALLER),
                                                       cycle_id,
                                                       curr->label().index(),
                                                       false)));
        }

        cycle.erase(curr);
        cycle.erase(next);
        curr = after_next;

        if (after_next != cycle.end()) {
          next = curr;
          ++next;
        }

        // Mark that we have reduced the cycle.
        reduced_cycle = true;
      }
      else {
        curr = next;
        if (next != cycle.end()) {
          ++next;
        }
      }

    } while (! is_last);

    // In case we have reduced the cycle, re-number the segments in it.
    if (reduced_cycle) {
      seg_index = 0;
      for (curr = cycle.begin(); curr != cycle.end(); ++curr)
        curr->label().set_index(seg_index++);
    }
    cycle.back().label().set_flag(true);
    return;
  }
};

} //namespace CGAL

#endif
