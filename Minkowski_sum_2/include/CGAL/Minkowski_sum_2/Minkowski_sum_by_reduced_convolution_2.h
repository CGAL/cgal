// Copyright (c) 2015  Tel-Aviv University (Israel).
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
// Author(s): Sebastian Morr    <sebastian@morr.cc>

#ifndef CGAL_MINKOWSKI_SUM_BY_REDUCED_CONVOLUTION_2_H
#define CGAL_MINKOWSKI_SUM_BY_REDUCED_CONVOLUTION_2_H

#include <CGAL/basic.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Arr_segment_traits_2.h>

#include <CGAL/Minkowski_sum_2/AABB_collision_detector_2.h>

#include <queue>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

namespace CGAL {

// This algorithm was first described by Evan Behar and Jyh-Ming Lien in "Fast
// and Robust 2D Minkowski Sum Using Reduced Convolution", IROS 2011.
// This implementation is based on Alon Baram's 2013 master's thesis "Polygonal
// Minkowski Sums via Convolution: Theory and Practice" at Tel-Aviv University.
template <class Kernel_, class Container_>
class Minkowski_sum_by_reduced_convolution_2
{
private:
  typedef Kernel_ Kernel;
  typedef Container_ Container;

  // Basic types:
  typedef CGAL::Polygon_2<Kernel, Container> Polygon_2;
  typedef CGAL::Polygon_with_holes_2<Kernel, Container> Polygon_with_holes_2;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Vector_2                     Vector_2;
  typedef typename Kernel::Direction_2                  Direction_2;
  typedef typename Kernel::Triangle_2                   Triangle_2;
  typedef typename Kernel::FT                           FT;

  // Segment-related types:
  typedef Arr_segment_traits_2<Kernel>                  Traits_2;
  typedef typename Traits_2::X_monotone_curve_2         Segment_2;
  typedef std::list<Segment_2>                          Segment_list;
  typedef Arr_default_dcel<Traits_2>                    Dcel;
  typedef std::pair<int, int>                           State;

  // Arrangement-related types:
  typedef Arrangement_with_history_2<Traits_2, Dcel>    Arrangement_history_2;
  typedef typename Arrangement_history_2::Halfedge_handle Halfedge_handle;
  typedef typename Arrangement_history_2::Face_iterator Face_iterator;
  typedef typename Arrangement_history_2::Face_handle   Face_handle;
  typedef typename Arrangement_history_2::Ccb_halfedge_circulator
  Ccb_halfedge_circulator;
  typedef typename Arrangement_history_2::Originating_curve_iterator
  Originating_curve_iterator;
  typedef typename Arrangement_history_2::Inner_ccb_iterator Inner_ccb_iterator;

  // Function object types:
  typename Kernel::Construct_translated_point_2 f_add;
  typename Kernel::Construct_vector_2 f_vector;
  typename Kernel::Construct_direction_2 f_direction;
  typename Kernel::Orientation_2 f_orientation;
  typename Kernel::Compare_xy_2 f_compare_xy;
  typename Kernel::Counterclockwise_in_between_2 f_ccw_in_between;

public:
  Minkowski_sum_by_reduced_convolution_2()
  {
    // Obtain kernel functors
    Kernel ker;
    f_add = ker.construct_translated_point_2_object();
    f_vector = ker.construct_vector_2_object();
    f_direction = ker.construct_direction_2_object();
    f_orientation = ker.orientation_2_object();
    f_compare_xy = ker.compare_xy_2_object();
    f_ccw_in_between = ker.counterclockwise_in_between_2_object();
  }

  template <typename OutputIterator>
  void operator()(const Polygon_2& pgn1, const Polygon_2& pgn2,
                  Polygon_2& outer_boundary, OutputIterator holes) const
  {
    CGAL_precondition(pgn1.is_simple());
    CGAL_precondition(pgn2.is_simple());
    CGAL_precondition(pgn1.orientation() == COUNTERCLOCKWISE);
    CGAL_precondition(pgn2.orientation() == COUNTERCLOCKWISE);

    const Polygon_with_holes_2 pwh1(pgn1);
    const Polygon_with_holes_2 pwh2(pgn2);

    common_operator(pwh1, pwh2, outer_boundary, holes);
  }

  template <class OutputIterator>
  void operator()(const Polygon_with_holes_2& pgn1,
                  const Polygon_with_holes_2& pgn2,
                  Polygon_2& outer_boundary, OutputIterator holes) const
  {
    common_operator(pgn1, pgn2, outer_boundary, holes);
  }

  template <class OutputIterator>
  void operator()(const Polygon_2& pgn1,
                  const Polygon_with_holes_2& pgn2,
                  Polygon_2& outer_boundary, OutputIterator holes) const
  {
    CGAL_precondition(pgn1.is_simple());
    CGAL_precondition(pgn1.orientation() == COUNTERCLOCKWISE);
    const Polygon_with_holes_2 pwh1(pgn1);
    common_operator(pwh1, pgn2, outer_boundary, holes);
  }

private:
  template <class OutputIterator>
  void common_operator(const Polygon_with_holes_2& pgn1,
                       const Polygon_with_holes_2& pgn2,
                       Polygon_2& outer_boundary, OutputIterator holes) const
  {
    // Initialize collision detector. It operates on pgn2 and on the inversed
    // pgn1:
    const Polygon_with_holes_2 inversed_pgn1 =
      transform(Aff_transformation_2<Kernel>(SCALING, -1), pgn1);
    AABB_collision_detector_2<Kernel, Container>
      collision_detector(pgn2, inversed_pgn1);

    // Compute the reduced convolution (see section 4.1 of Alon's master's
    // thesis)
    Segment_list reduced_convolution;
    build_reduced_convolution(pgn1, pgn2, reduced_convolution);

    // Insert the segments into an arrangement
    Arrangement_history_2 arr;
    insert(arr, reduced_convolution.begin(), reduced_convolution.end());

    // Trace the outer loop and put it in 'outer_boundary'
    get_outer_loop(arr, outer_boundary);

    // Check for each face whether it is a hole in the M-sum. If it is, add it
    // to 'holes'. See chapter 3 of of Alon's master's thesis.
    for (Face_iterator face = arr.faces_begin(); face != arr.faces_end();
         ++face)
    {
      handle_face(arr, face, holes, collision_detector);
    }
  }

  // Builds the reduced convolution for each pair of loop in the two
  // polygons-with-holes.
  void build_reduced_convolution(const Polygon_with_holes_2& pgnwh1,
                                 const Polygon_with_holes_2& pgnwh2,
                                 Segment_list& reduced_convolution) const
  {
    for (std::size_t x = 0; x < 1+pgnwh1.number_of_holes(); ++x)
    {
      for (std::size_t y = 0; y < 1+pgnwh2.number_of_holes(); ++y)
      {
        if ((x != 0) && (y != 0))
        {
          continue;
        }

        Polygon_2 pgn1, pgn2;

        if (x == 0) {
          pgn1 = pgnwh1.outer_boundary();
        }
        else {
          typename Polygon_with_holes_2::Hole_const_iterator it1 =
            pgnwh1.holes_begin();
          for (std::size_t count = 0; count < x-1; count++) { it1++; }
          pgn1 = *it1;
        }

        if (y == 0) {
          pgn2 = pgnwh2.outer_boundary();
        }
        else {
          typename Polygon_with_holes_2::Hole_const_iterator it2 =
            pgnwh2.holes_begin();
          for (std::size_t count = 0; count < y-1; count++) { it2++; }
          pgn2 = *it2;
        }

        build_reduced_convolution(pgn1, pgn2, reduced_convolution);
      }
    }
  }

  // Builds the reduced convolution using a fiber grid approach. For each
  // starting vertex, try to add two outgoing next states. If a visited
  // vertex is reached, then do not explore further. This is a BFS-like
  // iteration beginning from each vertex in the first column of the fiber
  // grid.
  void build_reduced_convolution(const Polygon_2& pgn1, const Polygon_2& pgn2,
                                 Segment_list& reduced_convolution) const
  {
    int n1 = static_cast<int>(pgn1.size());
    int n2 = static_cast<int>(pgn2.size());

    std::vector<Point_2> p1_vertices = vertices_of_polygon(pgn1);
    std::vector<Point_2> p2_vertices = vertices_of_polygon(pgn2);

    // Init the direcions of both polygons
    std::vector<Direction_2> p1_dirs = directions_of_polygon(p1_vertices);
    std::vector<Direction_2> p2_dirs = directions_of_polygon(p2_vertices);

    // Contains states that were already visited
    boost::unordered_set<State> visited_states;

    // Init the queue with vertices from the first column
    std::queue<State> state_queue;
    for (int i = n1-1; i >= 0; --i)
    {
      state_queue.push(State(i, 0));
    }

    while (state_queue.size() > 0)
    {
      State curr_state = state_queue.front();
      state_queue.pop();

      int i1 = curr_state.first;
      int i2 = curr_state.second;

      // If this state was already visited, skip it
      if (visited_states.count(curr_state) > 0)
      {
        continue;
      }
      visited_states.insert(curr_state);

      int next_i1 = (i1+1) % n1;
      int next_i2 = (i2+1) % n2;
      int prev_i1 = (n1+i1-1) % n1;
      int prev_i2 = (n2+i2-1) % n2;

      // Try two transitions: From (i,j) to (i+1,j) and to (i,j+1). Add
      // the respective segments, if they are in the reduced convolution.
      for(int step_in_pgn1 = 0; step_in_pgn1 <= 1; step_in_pgn1++)
      {
        int new_i1, new_i2;
        if (step_in_pgn1)
        {
          new_i1 = next_i1;
          new_i2 = i2;
        }
        else
        {
          new_i1 = i1;
          new_i2 = next_i2;
        }

        // If the segment's direction lies counterclockwise in between
        // the other polygon's vertex' ingoing and outgoing directions,
        // the segment belongs to the full convolution.
        bool belongs_to_convolution;
        if (step_in_pgn1)
        {
          belongs_to_convolution =
            f_ccw_in_between(p1_dirs[i1], p2_dirs[prev_i2], p2_dirs[i2]) ||
            p1_dirs[i1] == p2_dirs[i2];
        }
        else
        {
          belongs_to_convolution =
            f_ccw_in_between(p2_dirs[i2], p1_dirs[prev_i1], p1_dirs[i1]) ||
            p2_dirs[i2] == p1_dirs[prev_i1];
        }

        if (belongs_to_convolution)
        {
          state_queue.push(State(new_i1, new_i2));

          // Only edges added to convex vertices can be on the M-sum's boundary.
          // This filter only leaves the *reduced* convolution.
          bool convex;
          if (step_in_pgn1)
          {
            convex = is_convex(p2_vertices[prev_i2], p2_vertices[i2],
              p2_vertices[next_i2]);
          }
          else
          {
            convex = is_convex(p1_vertices[prev_i1], p1_vertices[i1],
              p1_vertices[next_i1]);
          }

          if (convex)
          {
            Point_2 start_point = get_point(i1, i2, p1_vertices, p2_vertices);
            Point_2 end_point = get_point(new_i1, new_i2, p1_vertices,
                                          p2_vertices);
            reduced_convolution.push_back(Segment_2(start_point, end_point));
          }
        }
      }
    }
  }

  // Returns a vector of the polygon's vertices, in case that Container
  // is std::list and we cannot use vertex(i).
  std::vector<Point_2> vertices_of_polygon(const Polygon_2& p) const
  {
    std::vector<Point_2> vertices;

    for (typename Polygon_2::Vertex_const_iterator it = p.vertices_begin();
         it != p.vertices_end(); it++)
    {
      vertices.push_back(*it);
    }
    return vertices;
  }

  // Returns a sorted list of the polygon's edges
  std::vector<Direction_2> directions_of_polygon(
    const std::vector<Point_2>& points) const
  {
    std::vector<Direction_2> directions;
    std::size_t n = points.size();

    for (std::size_t i = 0; i < n-1; ++i)
    {
      directions.push_back(f_direction(f_vector(points[i], points[i+1])));
    }
    directions.push_back(f_direction(f_vector(points[n-1], points[0])));

    return directions;
  }

  bool is_convex(const Point_2& prev, const Point_2& curr,
                 const Point_2& next) const
  {
    return f_orientation(prev, curr, next) == LEFT_TURN;
  }

  // Returns the point corresponding to a state (i,j).
  Point_2 get_point(int i1, int i2, const std::vector<Point_2>& pgn1,
                    const std::vector<Point_2>& pgn2) const
  {

    return f_add(pgn1[i1], Vector_2(Point_2(ORIGIN), pgn2[i2]));
  }

  // Put the outer loop of the arrangement in 'outer_boundary'
  void get_outer_loop(Arrangement_history_2& arr,
                      Polygon_2& outer_boundary) const
  {
    Inner_ccb_iterator icit = arr.unbounded_face()->inner_ccbs_begin();
    Ccb_halfedge_circulator circ_start = *icit;
    Ccb_halfedge_circulator circ = circ_start;

    do
    {
      outer_boundary.push_back(circ->source()->point());
    }
    while (--circ != circ_start);
  }

  // Check whether the face is on the M-sum's border. Add it to 'holes' if it is.
  template <class OutputIterator>
  void handle_face(const Arrangement_history_2& arr, const Face_handle face,
                   OutputIterator holes,
                   AABB_collision_detector_2<Kernel, Container>&
                   collision_detector) const
  {

    // If the face contains holes, it can't be on the Minkowski sum's border
    if (face->holes_begin() != face->holes_end())
    {
      return;
    }

    Ccb_halfedge_circulator start = face->outer_ccb();
    Ccb_halfedge_circulator circ = start;

    // The face needs to be orientable
    do
    {
      if (!do_original_edges_have_same_direction(arr, circ))
      {
        return;
      }
    }
    while (++circ != start);

    // When the reversed polygon 1, translated by a point inside of this face,
    // collides with polygon 2, this cannot be a hole
    Point_2 inner_point = get_point_in_face(face);
    if (collision_detector.check_collision(inner_point))
    {
      return;
    }

    // At this point, the face is a real hole, add it to 'holes'
    Polygon_2 pgn_hole;
    circ = start;

    do
    {
      pgn_hole.push_back(circ->source()->point());
    }
    while (--circ != start);

    *holes = pgn_hole;
    ++holes;
  }

  // Check whether the convolution's original edge(s) had the same direction as
  // the arrangement's half edge
  bool do_original_edges_have_same_direction(const Arrangement_history_2& arr,
                                             const Halfedge_handle he) const
  {
    Originating_curve_iterator segment_itr;

    for (segment_itr = arr.originating_curves_begin(he);
         segment_itr != arr.originating_curves_end(he); ++segment_itr)
    {
      if (f_compare_xy(segment_itr->source(), segment_itr->target()) ==
          (Comparison_result)he->direction())
      {
        return false;
      }
    }

    return true;
  }

  // Return a point in the face's interior by finding a diagonal
  Point_2 get_point_in_face(const Face_handle face) const
  {
    Ccb_halfedge_circulator current_edge = face->outer_ccb();
    Ccb_halfedge_circulator next_edge = current_edge;
    next_edge++;

    Point_2 a, v, b;

    // Move over the face's vertices until a convex corner is encountered:
    do
    {
      a = current_edge->source()->point();
      v = current_edge->target()->point();
      b = next_edge->target()->point();

      current_edge++;
      next_edge++;
    }
    while (!is_convex(a, v, b));

    Triangle_2 ear(a, v, b);
    FT min_distance = -1;
    const Point_2* min_q = 0;

    // Of the remaining vertices, find the one inside of the "ear" with minimal
    // distance to v:
    while (++next_edge != current_edge)
    {
      const Point_2& q = next_edge->target()->point();
      if (ear.has_on_bounded_side(q))
      {
        FT distance = squared_distance(q, v);
        if ((min_q == 0) || (distance < min_distance))
        {
          min_distance = distance;
          min_q = &q;
        }
      }
    }

    // If there was no vertex inside of the ear, return it's centroid.
    // Otherwise, return a point between v and min_q.
    if (min_q == 0)
    {
      return centroid(ear);
    }
    else
    {
      return midpoint(v, *min_q);
    }
  }

  template <class Transformation>
  Polygon_with_holes_2 transform(const Transformation& t,
                                 const Polygon_with_holes_2& p) const
  {
    Polygon_with_holes_2 result(CGAL::transform(t, p.outer_boundary()));

    typename Polygon_with_holes_2::Hole_const_iterator it = p.holes_begin();
    while(it != p.holes_end())
    {
      Polygon_2 p2(it->vertices_begin(), it->vertices_end());
      result.add_hole(CGAL::transform(t, p2));
      ++it;
    }
    return result;
  }
};

} // namespace CGAL

#endif
