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
//             (based on an old version by Eyal Flato)
//             Efi Fogel  <efifogel@gmail.com>

#ifndef CGAL_SMALL_SIDE_ANGLE_BISECTOR_DECOMPOSITION_2_H
#define CGAL_SMALL_SIDE_ANGLE_BISECTOR_DECOMPOSITION_2_H

#include <CGAL/Polygon_2.h>
#include <vector>
#include <list>

namespace CGAL {

/*!
 * \class
 * Small-side angle-bisector decomposition strategy.
 */
template <typename Kernel_,
          typename Container_ = std::vector<typename Kernel_::Point_2> >
class Small_side_angle_bisector_decomposition_2 {
public:
  typedef Kernel_                                        Kernel;
  typedef CGAL::Polygon_2<Kernel, Container_>            Polygon_2;
  typedef typename Kernel::Point_2                       Point_2;

private:
  typedef typename Kernel::Direction_2                   Direction_2;
  typedef typename Kernel::Segment_2                     Segment_2;

  typedef typename Polygon_2::Vertex_circulator          Vertex_circulator;

  // Kernel functors:
  typedef typename Kernel::Equal_2                       Equal_2;
  typedef typename Kernel::Construct_translated_point_2  Translate_point_2;
  typedef typename Kernel::Construct_vector_2            Construct_vector_2;
  typedef typename Kernel::Construct_direction_2         Construct_direction_2;
  typedef typename Kernel::Orientation_2                 Compute_orientation_2;
  typedef typename Kernel::Do_intersect_2                Do_intersect_2;
  typedef typename Kernel::Counterclockwise_in_between_2 Ccw_in_between_2;

  typedef std::list<unsigned int>                        Indices_list;
  typedef typename Indices_list::const_iterator          Indices_iterator;

  // Point with additional infomaration.
  struct Point_info_2 {
    Point_2           point;
    bool              is_reflex;
    unsigned int      reflex_count;
    Indices_list      visible;
    Indices_list      non_visible;

    // Default constructor.
    Point_info_2() :
      is_reflex(false),
      reflex_count(0),
      visible(),
      non_visible()
    {}

    // Check if the given vertex is visible.
    bool is_visible(unsigned int index) const
    {
      Indices_iterator it;
      for (it = visible.begin(); it != visible.end(); ++it) {
        if (*it == index) return (true);
      }
      return (false);
    }

    // Check if the given vertex is non-visible.
    bool is_non_visible(unsigned int index) const
    {
      Indices_iterator it;
      for (it = non_visible.begin(); it != non_visible.end(); ++it) {
        if (*it == index) return (true);
      }
      return (false);
    }
  };

  typedef std::vector<Point_info_2>                Point_vector_2;

  const Kernel* m_kernel;
  bool m_own_kernel;    // inidicates whether the kernel should be freed up.

  // Data members:
  Equal_2                 f_equal;
  Translate_point_2       f_add;
  Construct_vector_2      f_vector;
  Construct_direction_2   f_direction;
  Compute_orientation_2   f_orientation;
  Do_intersect_2          f_do_intersect;
  Ccw_in_between_2        f_ccw_in_between;

public:
  /*! Default constructor. */
  Small_side_angle_bisector_decomposition_2() :
    m_kernel(new Kernel),
    m_own_kernel(true)
  { init(); }

  /*! Constructor. */
  Small_side_angle_bisector_decomposition_2(const Kernel& kernel) :
    m_kernel(&kernel),
    m_own_kernel(false)
  { init(); }

  /*! Destructor */
  ~Small_side_angle_bisector_decomposition_2()
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
    f_orientation = m_kernel->orientation_2_object();
    f_do_intersect = m_kernel->do_intersect_2_object();
    f_ccw_in_between = m_kernel->counterclockwise_in_between_2_object();
  }

  /*!
   * Decompose a simple polygon to convex sub-polygons.
   * \param pgn The input polygon.
   * \param oi An output iterator of convex polygons.
   * \return A past-the-end iterator for the sub-polygons.
   */
  template <typename OutputIterator>
  OutputIterator operator()(const Polygon_2& pgn, OutputIterator oi) const
  {
    // Construct a point-info vector that represents the input polygon.
    Point_vector_2   vec;

    if (_construct_point_vector(pgn, vec)) {
      // The input polygon is convex - just return it as is (we use the
      // auxiliary function to make sure the returned polygon is oriented
      // counterclockwise).
      oi = _output_polygon(vec, oi);
      return (oi);
    }

    // Prepare a queue of polygons, represented as point vectors, to decompose.
    // We first try the small-side decomposition, and if we fail we turn to
    // the normal angle-bisector decomposition.
    std::list<Point_vector_2> ss_queue;
    std::list<Point_vector_2> ab_queue;
    unsigned int ind1, ind2;
    Point_vector_2 vec1, vec2;

    ind1 = ind2 = static_cast<unsigned int>(-1);        // Pacify some compiler

    ss_queue.push_back(vec);

    // Use the small-side decomposition as long as possible.
    while (! ss_queue.empty()) {
      // Extract the first polygon from the queue.
      vec = ss_queue.front();
      ss_queue.pop_front();

      // Try computing a minimal diagonal eliminating two reflex vertices.
      if (_compute_min_diagonal(vec, ind1, ind2)) {
        // We have found a minimal diagonal - use it to split the polygon.
        _split_by_diagonal(vec, ind1, ind2, false, vec1, vec2);

        // For either of the resulting sub-polygons: Report them if they are
        // convex, otherwise decompose them recursively.
        if (vec1[vec1.size() - 1].reflex_count == 0)
          oi = _output_polygon(vec1, oi);
        else
          ss_queue.push_back(vec1);

        if (vec2[vec2.size() - 1].reflex_count == 0)
          oi = _output_polygon(vec2, oi);
        else
          ss_queue.push_back(vec2);
      }
      else {
        // We have failed to find a minimal diagonal eliminating two reflex
        // vertices - use the simple angle bisector method for this polygon.
        ab_queue.push_back(vec);
      }
    }

    // Decompose the remaining polygons.
    while (! ab_queue.empty()) {
      // Extract the first polygon from the queue.
      vec = ab_queue.front();
      ab_queue.pop_front();

      // Compute the best angle bisector and split the polygon accordingly.
      _approximate_angle_bisector(vec, ind1, ind2);

      _split_by_diagonal(vec, ind1, ind2, true, vec1, vec2);

      // For either of the resulting sub-polygons: Report them if they are
      // convex, otherwise decompose them recursively.
      if (vec1[vec1.size() - 1].reflex_count == 0)
        oi = _output_polygon(vec1, oi);
      else
        ab_queue.push_back(vec1);

      if (vec2[vec2.size() - 1].reflex_count == 0)
        oi = _output_polygon(vec2, oi);
      else
        ab_queue.push_back(vec2);
    }

    return (oi);
  }

private:

  /*! Return the succesive index of a 'point info' vector. */
  inline unsigned int _vec_succ(const Point_vector_2& vec,
                                unsigned int i) const
  {
    return ((i + 1) % vec.size());
  }

  /*! Return the previous index of a 'point info' vector. */
  inline unsigned int _vec_pred(const Point_vector_2& vec,
                                unsigned int i) const
  {
    if (i == 0) return static_cast<unsigned int>(vec.size() - 1);
    return (i - 1);
  }

  /*!
   * Reflect the point b around the point a.
   * \return The result of: a + vec('a'-'b').
   */
  Point_2 _opposite_from_vertex(const Point_2& a, const Point_2& b) const
  {
    return (f_add(a, f_vector(b, a)));
  }

  // The "reflex zone" of a vertex v is the slice created by continuing the
  // edges incident to this vertex, from v to infinity. This enumeration
  // describes the location of a point with respect to this zone.
  enum Reflex_zone_position {
    LEFT_OF_ZONE = -1,
    INSIDE_ZONE = 0,
    RIGHT_OF_ZONE = 1,
    ON_ZONE_BOUNDARY = 3,
    OPPOSITE_OF_ZONE = 9
  };

  /*!
   * Get the position of a query point with respect to the reflex zone of a
   * given vertex.
   * \param vec A vector defining counterclockwise-oriented polygon.
   * \param v_ind The index of the vertex v.
   * \param q_pt A query point.
   * \return The location of the query point with respect to the zone.
   */
  Reflex_zone_position _query_reflex_zone(const Point_vector_2& vec,
                                          unsigned int v_ind,
                                          const Point_2& q_pt) const
  {
    // Initialize q, the reflex vertex, and it previous and next vertices
    // and incident edges.
    const Point_2&  v_pt = vec[v_ind].point;
    const Point_2&  pred_pt = vec[_vec_pred(vec, v_ind)].point;
    const Point_2&  succ_pt = vec[_vec_succ(vec, v_ind)].point;
    Direction_2     pred_dir = f_direction(f_vector(v_pt, pred_pt));
    Direction_2     succ_dir = f_direction(f_vector(v_pt, succ_pt));
    Point_2         opp_pred_pt = _opposite_from_vertex(v_pt, pred_pt);
    Point_2         opp_succ_pt = _opposite_from_vertex(v_pt, succ_pt);
    Direction_2     opp_pred_dir = f_direction(f_vector(v_pt, opp_pred_pt));
    Direction_2     opp_succ_dir = f_direction(f_vector(v_pt, opp_succ_pt));
    Direction_2     q_dir = f_direction(f_vector(v_pt, q_pt));

    if (f_ccw_in_between(q_dir, opp_pred_dir, opp_succ_dir)) {
      // The point is contained between the continuation of the edges around v,
      // so it is inside the reflex zone.
      return (INSIDE_ZONE);
    }

    if (f_equal(q_dir, opp_pred_dir) || f_equal(q_dir, opp_succ_dir)) {
      // The point lies on the boundary of the reflex zone.
      return (ON_ZONE_BOUNDARY);
    }

    if (f_ccw_in_between(q_dir, opp_succ_dir, pred_dir)) {
      // The point lies to the left of the reflex zone.
      return (LEFT_OF_ZONE);
    }

    if (f_ccw_in_between(q_dir, succ_dir, opp_pred_dir)) {
      // The point lies to the right of the reflex zone.
      return (RIGHT_OF_ZONE);
    }

    // If we reached here, the query point lies on the opposite side of the
    // reflex zone.
    CGAL_assertion(f_ccw_in_between(q_dir, pred_dir, succ_dir) ||
                   f_equal(q_dir, pred_dir) || f_equal(q_dir, succ_dir));

    return (OPPOSITE_OF_ZONE);
  }

  /*!
   * Check whether a given vertex u is visible from a reflex vertex v.
   * \param vec A vector defining counterclockwise-oriented polygon.
   * \param v_ind The index of the reflex vertex v.
   * \param u_ind The index of the vertex u.
   * \return Whether the vertices are visible.
   */
  bool _is_visible(Point_vector_2& vec,
                   unsigned int v_ind, unsigned int u_ind) const
  {
    CGAL_precondition(vec[v_ind].is_reflex);

    // Check whether the visiblity status is already known.
    if (vec[v_ind].is_visible(u_ind)) return (true);
    if (vec[v_ind].is_non_visible(u_ind)) return (false);

    // Check whether the diagonal (v_ind, u_ind) does not intersect any of
    // the polygn edges.
    bool    result;

    if ((u_ind == v_ind) ||
        (u_ind == _vec_pred(vec, v_ind)) || (u_ind == _vec_succ(vec, v_ind)))
    {
      // A polygon edge (or a degenerate diagonal) is not considered as a legal
      // visibility diagonal.
      result = false;
    }
    else if (_query_reflex_zone(vec, v_ind,
                                vec[u_ind].point) == OPPOSITE_OF_ZONE)
    {
      // In this case the diagonal between u and v lies outside the polygon,
      // so the two vertices are not visible.
      result = false;
    }
    else {
      // Perform an exhaustive search on the polygon edges and check if any of
      // them intersects the diagonal segment uv.
      Segment_2 diag(vec[u_ind].point, vec[v_ind].point);
      unsigned int curr = _vec_succ(vec, v_ind);
      unsigned int next = _vec_succ(vec, curr);

      result = true;
      while (result && next != v_ind) {
        if (curr != u_ind && next != u_ind) {
          if (do_intersect(diag, Segment_2(vec[curr].point, vec[next].point)))
            result = false;
        }

        curr = next;
        next = _vec_succ(vec, curr);
      }
    }

    // Update the visibility status.
    if (result) vec[v_ind].visible.push_back(u_ind);
    else vec[v_ind].non_visible.push_back(u_ind);

    return (result);
  }

  /*!
   * Construct a Point_info_2 vector from the given polygon.
   * \param poly The input polygon.
   * \param vec Output: The point vector.
   * \return Whether the polygon is convex.
   */
  bool _construct_point_vector(const Polygon_2& pgn, Point_vector_2& vec) const
  {
    // Resize the vector to fit the polygon.
    const unsigned int n = static_cast<unsigned int>(pgn.size());
    const bool forward = (pgn.orientation() == COUNTERCLOCKWISE);
    Vertex_circulator prev, curr, next;
    unsigned int reflex_count = 0;
    unsigned int k;

    vec.resize(n);

    // Initialize the vertex circulators.
    next = curr = prev = pgn.vertices_circulator();

    if (forward) --prev;
    else ++prev;

    // Traverse the polygon's vertices in a counterclockwise order and prepare
    // the output vector.
    for (k = 0; k < n; ++k) {
      // Set the next vertex.
      if (forward) ++next;
      else --next;

      // Set the current point.
      vec[k].point = *curr;

      // Check if the current vertex is reflex by checking the orientation of
      // the polygon around it.
      vec[k].is_reflex = (f_orientation(*prev, *curr, *next) == RIGHT_TURN);

      // Set the number of reflex vertices from the beginning of the vector
      // until the k'th (including itself).
      if (vec[k].is_reflex) ++reflex_count;

      vec[k].reflex_count = reflex_count;

      // Proceed to the next vertex.
      prev = curr;
      curr = next;
    }

    // The polygon is convex if there are not reflex vertices:
    return (reflex_count == 0);
  }

  /*!
   * Given a polygon and to of its vertices u and v, count the number of reflex
   * vertices between u and v, and between v and u (not including u and v
   * themselves), and return the lesser quantity.
   * \param vec A vector defining counterclockwise-oriented polygon.
   * \param v_ind The index of the vertex v.
   * \param u_ind The index of the vertex u.
   * \return The minimal number of reflex vertices between u and v.
   */
  unsigned int _count_reflex_vertices(const Point_vector_2& vec,
                                      unsigned int u_ind,
                                      unsigned int v_ind) const
  {
    if (u_ind == v_ind) return (0);

    unsigned int count1 = 0, count2 = 0;
    unsigned int k;

    // Make sure that v_ind > u_ind (if necessary, swap indices).
    if (u_ind > v_ind) {
      k = u_ind;
      u_ind = v_ind;
      v_ind = k;
    }

    // Compute how many reflex vertices there are between u and v (not
    // including v itself).
    count1 = (vec[v_ind].reflex_count - vec[u_ind].reflex_count);

    if (vec[v_ind].is_reflex) --count1;

    // Compute how many reflex vertices there are between v and u (not
    // including u itself).
    count2 = (vec[vec.size() - 1].reflex_count - vec[v_ind].reflex_count) +
      vec[u_ind].reflex_count;

    if (vec[u_ind].is_reflex) --count2;

    // Return the smaller value.
    if (count1 < count2) return count1;
    else return count2;
  }

  /*!
   * Given a non-convex polygon, try to find a diagonal connecting to mutually
   * visible reflex vertices u and v such that the number of reflex vertices
   * between u and v is minimized.
   * \param vec A vector defining counterclockwise-oriented polygon.
   * \param v_ind Output: The index of the vertex v.
   * \param u_ind Output: The index of the vertex u.
   * \return Whether a diagonal was found.
   */
  bool _compute_min_diagonal(Point_vector_2& vec,
                             unsigned int& u_ind, unsigned int& v_ind) const
  {
    unsigned int ind1, ind2;
    Reflex_zone_position zone_pos;
    unsigned int curr_count;
    unsigned int min_count = 0;
    bool found = false;
    unsigned int dist;

    // Let the distance between pairs of vertices we try vary between 2 and
    // half the size of the polygon.
    for (dist = 2; dist < (vec.size() + 1) / 2; ++dist) {
      // Traverse the polygon vertices.
      for (ind1 = 0; ind1 < vec.size(); ++ind1) {
        // Ignore convex vertices.
        if (! vec[ind1].is_reflex) continue;

        ind2 = (ind1 + dist) % vec.size();

        // Ignore convex vertices.
        if (! vec[ind2].is_reflex) continue;

        // Count the number of reflex vertices between the current pair of
        // reflex vertices.
        curr_count = _count_reflex_vertices(vec, ind1, ind2);

        // If we are not below the minimal, stop checking.
        if (found && (min_count < curr_count)) continue;

        // Check whether the two vertices lie in the reflex zone of one
        // another. If so, the diagonal between them splits each of the angles
        // of u and v into two angles than are not greater than 180 degrees.
        zone_pos = _query_reflex_zone(vec, ind1, vec[ind2].point);

        if ((zone_pos != INSIDE_ZONE) && (zone_pos != ON_ZONE_BOUNDARY))
          continue;

        zone_pos = _query_reflex_zone(vec, ind2, vec[ind1].point);

        if ((zone_pos != INSIDE_ZONE) && (zone_pos != ON_ZONE_BOUNDARY))
          continue;

        // Check whether the two vertices are mutually visible.
        if (! _is_visible(vec, ind1, ind2)) continue;

        // If we passed all the various tests and reached here, we can update
        // the minimal count.
        min_count = curr_count;
        u_ind = ind1;
        v_ind = ind2;
        found = true;
      }
    }

    // Return whether a diagonal has been found.
    return (found);
  }

  /*!
   * Split the given polygon by the given diagonal.
   * \param vec A vector defining counterclockwise-oriented polygon.
   * \param ind1 The index of the first diagonal end-vertex.
   * \param ind2 The index of the second diagonal end-vertex.
   * \param check_ends Should we check for convexity near ind1 and ind2.
   * \param vec1 Output: The first output polygon.
   * \param vec2 Output: The second output polygon.
   */
  void _split_by_diagonal(const Point_vector_2& vec,
                          unsigned int ind1, unsigned int ind2,
                          bool check_ends,
                          Point_vector_2& vec1,
                          Point_vector_2& vec2) const
  {
    unsigned int i1, i2;
    Indices_iterator iter;
    unsigned int reflex_count;

    // Make sure that the index ind2 is greater than ind1.
    if (ind1 > ind2) {
      i1 = ind1;
      ind1 = ind2;
      ind2 = i1;
    }

    // Resize the output vectors.
    vec1.resize(ind2 - ind1 + 1);
    vec2.resize(vec.size() - (ind2 - ind1 - 1));

    // Construct the first output polygon. Notice we traverse the polygon in
    // a counterclockwise direction from ind1 to ind2 and back.
    reflex_count = 0;

    for (i1 = ind1; i1 <= ind2; ++i1) {
      // Copy the point.
      vec1[i1 - ind1].point = vec[i1].point;

      // Clear the visible vertices list.
      vec1[i1-ind1].visible.clear();
      vec1[i1-ind1].non_visible.clear();

      if (i1 == ind1) {
        if (check_ends) {
          // Check the turn around ind1.
          vec1[i1-ind1].is_reflex =
            (f_orientation(vec[ind2].point,
                            vec[ind1].point,
                            vec[_vec_succ(vec, ind1)].point) == RIGHT_TURN);
        }
        else {
          // The two diagonal end-vertices are not reflex any more.
          vec1[i1-ind1].is_reflex = false;
        }
      }
      else if (i1 == ind2) {
        if (check_ends) {
          // Check the turn around ind2.
          vec1[i1-ind1].is_reflex =
            (f_orientation(vec[_vec_pred(vec, ind2)].point,
                            vec[ind2].point,
                            vec[ind1].point) == RIGHT_TURN);
        }
        else {
          // The two diagonal end-vertices are not reflex any more.
          vec1[i1-ind1].is_reflex = false;
        }
      }
      else {
        // Copy the reflexity flag.
        vec1[i1-ind1].is_reflex = vec[i1].is_reflex;

        // Set the visible vertices list.
        for (iter = vec[i1].visible.begin();
             iter != vec[i1].visible.end(); ++iter)
        {
          if ((*iter > ind1) && (*iter < ind2))
            vec1[i1-ind1].visible.push_back(*iter - ind1);
        }

        // Set the non-visible vertices list.
        for (iter = vec[i1].non_visible.begin();
             iter != vec[i1].non_visible.end(); ++iter)
        {
          if ((*iter > ind1) && (*iter < ind2))
            vec1[i1-ind1].non_visible.push_back(*iter - ind1);
        }
      }

      // Set the reflex count.
      if (vec1[i1-ind1].is_reflex)
        reflex_count++;

      vec1[i1-ind1].reflex_count = reflex_count;
    }

    // Construct the second output polygon. Notice we traverse the polygon in
    // a counterclockwise direction from ind2 to ind1 and back.
    reflex_count = 0;
    i2 = 0;
    i1 = ind2;

    // go ccw on all the vertices from 'ind2' to 'ind1'.
    while (i1 != (ind1+1)) {
      // Copy the point.
      vec2[i2].point = vec[i1].point;

      // Clear the visible vertices list.
      vec2[i2].visible.clear();
      vec2[i2].non_visible.clear();

      if (i1 == ind1) {
        if (check_ends) {
          // Check the turn around ind1.
          vec2[i2].is_reflex =
            (f_orientation(vec[_vec_pred(vec, ind1)].point,
                            vec[ind1].point,
                            vec[ind2].point) == RIGHT_TURN);
        }
        else {
          // The two diagonal end-vertices are not reflex any more.
          vec2[i2].is_reflex = false;
        }
      }
      else if (i1 == ind2) {
        if (check_ends) {
          // Check the turn around ind2.
          vec2[i2].is_reflex =
            (f_orientation(vec[ind1].point,
                            vec[ind2].point,
                            vec[_vec_succ(vec, ind2)].point) == RIGHT_TURN);
        }
        else {
          // The two diagonal end-vertices are not reflex any more.
          vec2[i2].is_reflex = false;
        }
      }
      else {
        // Copy the reflexity flag.
        vec2[i2].is_reflex = vec[i1].is_reflex;

        // Set the visible vertices list.
        for (iter = vec[i1].visible.begin();
             iter != vec[i1].visible.end(); ++iter)
        {
          if (*iter > ind2) {
            vec2[i2].visible.push_back(*iter - ind2);
          }
          else if (*iter < ind1) {
            vec2[i2].visible.push_back(*iter +
                                       static_cast<unsigned int>(vec.size()) -
                                       ind2);
          }
        }

        // Set the non-visible vertices list.
        for (iter = vec[i1].non_visible.begin();
             iter != vec[i1].non_visible.end(); ++iter)
        {
          if (*iter > ind2) {
            vec2[i2].non_visible.push_back(*iter - ind2);
          }
          else if (*iter < ind1) {
            vec2[i2].non_visible.push_back(*iter +
                                           static_cast<int>(vec.size()) -
                                           ind2);
          }
        }
      }

      // Set the reflex count.
      if (vec2[i2].is_reflex)
        reflex_count++;

      vec2[i2].reflex_count = reflex_count;

      // advance the indices.
      i1 = _vec_succ(vec, i1);
      i2++;
    }

    return;
  }

  /*!
   * Get the angle ratio created by the the bisection of the angle at the
   * reflex vertex v by the diagonal uv.
   * \param vec A vector defining counterclockwise-oriented polygon.
   * \param v_ind The index of the vertex v.
   * \param u_ind The index of the vertex u.
   * \return The ratio between the bisected angles.
   */
  double _angle_bisection_ratio(Point_vector_2& vec,
                                unsigned int v_ind,
                                unsigned int u_ind) const
  {
    const double    _2_PI = 6.2831853;

    // Get the points at u and v, as well as v's adjacencies and approximate
    // their coordinates.
    const Point_2&  u_pt = vec[u_ind].point;
    const Point_2&  v_pt = vec[v_ind].point;
    const Point_2&  pred_pt = vec[_vec_pred(vec, v_ind)].point;
    const Point_2&  succ_pt = vec[_vec_succ(vec, v_ind)].point;
    const double    ux = CGAL::to_double(u_pt.x());
    const double    uy = CGAL::to_double(u_pt.y());
    const double    vx = CGAL::to_double(v_pt.x());
    const double    vy = CGAL::to_double(v_pt.y());
    const double    pred_x = CGAL::to_double(pred_pt.x());
    const double    pred_y = CGAL::to_double(pred_pt.y());
    const double    succ_x = CGAL::to_double(succ_pt.x());
    const double    succ_y = CGAL::to_double(succ_pt.y());

    // Compute the edge length and the diagonal length.
    const double len_uv =
      ::sqrt((ux - vx) * (ux - vx) + (uy - vy) * (uy - vy));
    const double len_pred =
      ::sqrt((pred_x - vx) * (pred_x - vx) + (pred_y - vy) * (pred_y - vy));
    const double len_succ =
      ::sqrt((succ_x - vx) * (succ_x - vx) + (succ_y - vy) * (succ_y - vy));

    // Compute the angle a1 = <) (pred_pt, v_pt, u_pt):
    const double cos_a1 = ((ux - vx) * (pred_x - vx) +
                           (uy - vy) * (pred_y - vy)) / (len_uv * len_pred);
    double a1 = ::acos(cos_a1);

    // The angle a1 is larger than 180 degree:
    if (f_orientation(pred_pt, v_pt, u_pt) == RIGHT_TURN) a1 = _2_PI - a1;

    // Compute the angle a2 = <) (u_pt, v_pt, succ_pt):
    const double cos_a2 = ((ux - vx) * (succ_x - vx) +
                           (uy - vy) * (succ_y - vy)) / (len_uv * len_succ);
    double a2 = ::acos(cos_a2);

    if (f_orientation(u_pt, v_pt, succ_pt) == RIGHT_TURN)
      // The angle a1 is larger than 180 degree:
      a2 = _2_PI - a2;

    // Return the ratio of max(a1,a2)/min(a1,a2).
    if (a1 == 0 || a2 == 0) return (10000);

    if (a1 > a2) return (a1 / a2);
    else return (a2 / a1);
  }

  /*!
   * Perform an approximate angle-bisector decomposition, by locating a reflex
   * vertex and another vertex such that the diagonal between them is not
   * blocked, and such that it approximates the angle bisector of the reflex
   * vertex.
   * \param vec A vector defining counterclockwise-oriented polygon.
   * \param reflex_ind Output: The index of the reflex vertex whose angle we
   *                           wish to bisect.
   * \param other_ind Output: The index of the other polygon vertex.
   */
  void _approximate_angle_bisector(Point_vector_2& vec,
                                   unsigned int& reflex_ind,
                                   unsigned int& other_ind) const
  {
    unsigned int ind1, ind2;
    double curr_ratio;
    double min_ratio = 0;
    bool found = false;
    unsigned int dist;

    // Let the distance between pairs of vertices we try vary between 2 and
    // the size of the polygon - 2.
    for (dist = 2; dist <= vec.size() - 2; ++dist) {
      // Traverse the polygon vertices.
      for (ind1 = 0; ind1 < vec.size(); ++ind1) {
        // Ignore convex vertices.
        if (! vec[ind1].is_reflex) continue;

        ind2 = (ind1 + dist) % vec.size();

        // Get the current ratio between the two angles the current
        // diagonal splits around vec[ind1].
        curr_ratio = _angle_bisection_ratio(vec, ind1, ind2);

        // If we are not below the minimal, stop checking.
        if (found && min_ratio < curr_ratio) continue;

        // Check whether the two vertices are mutually visible.
        if (! _is_visible(vec, ind1, ind2)) continue;

        // If we passed all the various tests and reached here, we can update
        // the minimal angle difference.
        min_ratio = curr_ratio;
        reflex_ind = ind1;
        other_ind = ind2;
        found = true;
      }
    }

    // Return whether a diagonal has been found.
    CGAL_assertion(found);
    return;
  }

  /*!
   * Convert the point-info vector into a polygon and send it as an output.
   * \param vec A vector defining counterclockwise-oriented convex polygon.
   * \param oi The output iterator.
   * \return A past-the-end iterator for the sub-polygons.
   */
  template <typename OutputIterator>
  OutputIterator _output_polygon(const Point_vector_2& vec,
                                 OutputIterator oi) const
  {
    const unsigned int n = static_cast<int>(vec.size());
    Polygon_2 pgn;
    for (unsigned int k = 0; k < n; ++k) pgn.push_back(vec[k].point);
    *oi++ = pgn;
    return (oi);
  }
};

} //namespace CGAL

#endif
