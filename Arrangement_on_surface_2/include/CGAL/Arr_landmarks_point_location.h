// Copyright (c) 2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
//
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
//                 Ron Wein     <wein@post.tau.ac.il>

#ifndef CGAL_ARR_LANDMARKS_POINT_LOCATION_H
#define CGAL_ARR_LANDMARKS_POINT_LOCATION_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * Definition of the Arr_landmarks_point_location<Arrangement> template.
 */

//#define CGAL_DEBUG_LM

#include <CGAL/Arr_point_location_result.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_point_location/Arr_lm_vertices_generator.h>

#include <set>

namespace CGAL {

/*! \class Arr_landmarks_point_location
 * A class that answers point-location queries on an arrangement using the
 * landmarks algorithm, namely by locating the (approximately) nearest
 * landmark point to the qury point and walking from it toward the query
 * point.
 * This class-template has two parameters:
 * Arrangement corresponds to an arrangement-on-surface instantiation.
 * Generator is a class that generates the set of landmarks.
 */

template <class Arrangement_,
          class Generator_ = Arr_landmarks_vertices_generator<Arrangement_> >
class Arr_landmarks_point_location
{
public:
  typedef Arrangement_                                  Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2     Geometry_traits_2;
  typedef Generator_                                    Generator;

  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

  typedef typename Arrangement_2::Vertex_const_iterator Vertex_const_iterator;
  typedef typename Arrangement_2::Halfedge_const_iterator
                                                        Halfedge_const_iterator;
  typedef typename Arrangement_2::Halfedge_around_vertex_const_circulator
    Halfedge_around_vertex_const_circulator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Outer_ccb_const_iterator
    Outer_ccb_const_iterator;
  typedef typename Arrangement_2::Inner_ccb_const_iterator
    Inner_ccb_const_iterator;
  typedef typename Arrangement_2::Isolated_vertex_const_iterator
    Isolated_vertex_const_iterator;

  typedef typename Arrangement_2::Point_2               Point_2;
  typedef typename Arrangement_2::X_monotone_curve_2    X_monotone_curve_2;

  typedef Arr_point_location_result<Arrangement_2>      Result;
  typedef typename Result::Type                         Result_type;

  // Support cpp11::result_of
  typedef Result_type                                   result_type;

protected:
  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2> Traits_adaptor_2;

  /*! \struct Less_halfedge_handle
   * Used to sort handles.
   */
  struct Less_halfedge_handle {
    bool operator()(Halfedge_const_handle h1, Halfedge_const_handle h2) const
    { return (&(*h1) < &(*h2)); }
  };

  typedef std::set<Halfedge_const_handle, Less_halfedge_handle> Halfedge_set;

  // Data members:
  const Arrangement_2*    p_arr;     // The associated arrangement.
  const Traits_adaptor_2* m_traits;  // Its associated traits object.
  Generator*              lm_gen;    // The associated landmark generator.
  bool                    own_gen;   // Indicates whether the generator
                                     // has been locally allocated.

  template<typename T>
  Result_type make_result(T t) const { return Result::make_result(t); }
  inline Result_type default_result() const { return Result::default_result(); }

public:
  /*! Default constructor. */
  Arr_landmarks_point_location() :
    p_arr(NULL),
    m_traits(NULL),
    lm_gen(NULL),
    own_gen(false)
  {}

  /*! Constructor given an arrangement only. */
  Arr_landmarks_point_location(const Arrangement_2& arr) :
    p_arr(&arr),
    m_traits(static_cast<const Traits_adaptor_2*>(p_arr->geometry_traits())),
    lm_gen(new Generator(arr)),         // allocate the landmarks generator.
    own_gen(true)
  { }

  /*! Constructor given an arrangement, and landmarks generator. */
  Arr_landmarks_point_location(const Arrangement_2& arr, Generator* gen) :
    p_arr(&arr),
    m_traits(static_cast<const Traits_adaptor_2*>(p_arr->geometry_traits())),
    lm_gen(gen),
    own_gen(false)
  { }

  /*! Destructor. */
  ~Arr_landmarks_point_location()
  {
    if (own_gen) {
      delete lm_gen;
      lm_gen = NULL;
    }
  }

 /*! Attach an arrangement object (and a generator, if supplied). */
  void attach(const Arrangement_2& arr, Generator* gen = NULL)
  {
    // Keep a pointer to the associated arrangement.
    p_arr = &arr;
    m_traits = static_cast<const Traits_adaptor_2*>(p_arr->geometry_traits());

    // Update the landmarks generator.
    if (gen != NULL) {
      // In case a generator is given, keep a pointer to it.
      CGAL_assertion(lm_gen == NULL);
      lm_gen = gen;
      own_gen = false;
    }
    else if (lm_gen != NULL) {
      // In case a generator exists internally, make sure it is attached to
      // the given arrangement.
      Arrangement_2& non_const_arr = const_cast<Arrangement_2&>(*p_arr);
      lm_gen->attach(non_const_arr);
    }
    else {
      // Allocate a new generator, attached to the given arrangement.
      lm_gen = new Generator(arr);
      own_gen = true;
    }
  }

  /*! Detach the instance from the arrangement object. */
  void detach()
  {
    p_arr = NULL;
    m_traits = NULL;

    CGAL_assertion(lm_gen != NULL);
    if (lm_gen)
      lm_gen->detach();
  }

  /*!
   * Locate the arrangement feature containing the given point.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type locate(const Point_2& p) const;

protected:
  /*! Walk from the given vertex to the query point.
   * \param vh The given vertex handle.
   * \param p The query point.
   * \param crossed_edges In/Out: The set of edges crossed so far.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type _walk_from_vertex(Vertex_const_handle vh,
                                const Point_2& p,
                                Halfedge_set& crossed_edges) const;

  /*! Locate an edge around a given vertex that is the predecessor of the
   * curve connecting the vertex to the query point in a clockwise order.
   * \param vh The vertex.
   * \param p The query point.
   * \param new_vertex Output: Whether a closer vertex to p was found.
   * \return The desired object (a halfedge handle or a vertex handle).
   */
  result_type _find_face_around_vertex(Vertex_const_handle vh,
                                       const Point_2& p,
                                       bool& new_vertex) const;

  /*! Walk from a point on a given halfedge to the query point.
   * \param eh The given halfedge handle.
   * \param np The point that the walk starts from.
   * \param p The query point.
   * \param crossed_edges In/Out: The set of edges crossed so far.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type _walk_from_edge(Halfedge_const_handle eh,
                              const Point_2& np,
                              const Point_2& p,
                              Halfedge_set& crossed_edges) const;
  /*! In case the arrangement's curve contained in the segment
   * from the nearest landmark to the query point
   * \param he The given halfedge handle.
   * \param p_is_left Is the query point the left endpoint of seg.
   * \param p The query point.
   * \param crossed_edges In/Out: The set of edges crossed so far.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type
  _deal_with_curve_contained_in_segment(Halfedge_const_handle he,
                                        bool p_is_left,
                                        const Point_2& p,
                                        Halfedge_set& crossed_edges) const;

  /*! Walk from a point in a face to the query point.
   * \param fh A halfedge handle that points to the face.
   * \param np The point that the walk starts from.
   * \param p The query point.
   * \param crossed_edges In/Out: The set of edges crossed so far.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type _walk_from_face(Face_const_handle fh,
                              const Point_2& np,
                              const Point_2& p,
                              Halfedge_set& crossed_edges) const;

  /*! Find a halfedge on the given CCB that intersects the given x-monotone
   * curve, connecting the current landmark to the query point.
   * \param circ The CCB circulator.
   * \param seg The segment connecting the landmark and the query point.
   * \param p The query point.
   * \param p_is_left Is the query point the left endpoint of seg.
   * \param crossed_edges In/Out: The set of edges crossed so far.
   * \param is_on_edge Output: Does the query point p lies on the edge.
   * \param is_target Output: Is the query point p equal to the target vertex.
   * \param new_vertex Output: if found a closer vertex to the query point.
   * \param cv_is_contained_in_seg Output: Whether cv is contained inside seg.
   * \return A handle to the halfedge (if no intersecting edge is found, the
   *         function returns an ivalid halfedge handle).
   */
  Halfedge_const_handle
  _intersection_with_ccb(Ccb_halfedge_const_circulator circ,
                         const X_monotone_curve_2& seg,
                         const Point_2& p,
                         bool p_is_left,
                         Halfedge_set& crossed_edges,
                         bool& is_on_edge,
                         bool& is_target,
                         bool& cv_is_contained_in_seg,
                         Vertex_const_handle& new_vertex) const;

  /*! Return the halfedge that contains the query point.
   * \param he The halfedge handle.
   * \param crossed_edges In/Out: The set of edges crossed so far.
   * \param p The query point.
   * \param is_target Output: Is the query point p equal to the target vertex.
   */
  Halfedge_const_handle
  _in_case_p_is_on_edge(Halfedge_const_handle he,
                        Halfedge_set& crossed_edges,
                        const Point_2& p,
                        bool& is_target) const;

  /*! Check whether the given curve intersects a simple segment, which connects
   * the current landmark to the query point, an odd number of times.
   * \param cv The curve.
   * \param seg The segment connecting the landmark and the query point.
   * \param p_is_left Is the query point the left endpoint of seg.
   * \param p_on_curve Output: Whether p lies on cv.
   * \param cv_and_seg_overlap Output: Whether cv and seg overlap.
   * \param cv_is_contained_in_seg Output: Whether cv is contained inside seg.
   * \return Whether the two curves have an odd number of intersections.
   */
   bool _have_odd_intersections(const X_monotone_curve_2& cv,
                                const X_monotone_curve_2& seg,
                                bool p_is_left,
                                bool& p_on_curve,
                                bool& cv_and_seg_overlap,
                                bool& cv_is_contained_in_seg) const;
};

} //namespace CGAL

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Arr_landmarks_pl_impl.h>

#include <CGAL/enable_warnings.h>

#endif
