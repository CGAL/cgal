// Copyright (c) 2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Idit Haran   <haranidi@post.tau.ac.il>
//             Ron Wein     <wein@post.tau.ac.il>

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
#include <CGAL/Arr_tags.h>

#include <set>

namespace CGAL {

/*! \class Arr_landmarks_point_location
 * A class that answers point-location queries on an arrangement using the
 * landmarks algorithm, namely by locating the (approximately) nearest
 * landmark point to the query point and walking from it toward the query
 * point.
 * This class-template has two parameters:
 * Arrangement corresponds to an arrangement-on-surface instantiation.
 * Generator is a class that generates the set of landmarks.
 */

template <typename Arrangement_,
          typename Generator_ = Arr_landmarks_vertices_generator<Arrangement_>>
class Arr_landmarks_point_location {
public:
  using Arrangement_2 = Arrangement_;
  using Generator = Generator_;
  using Geometry_traits_2 = typename Arrangement_2::Geometry_traits_2;

  using Vertex_const_handle = typename Arrangement_2::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement_2::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement_2::Face_const_handle;

  using Vertex_const_iterator = typename Arrangement_2::Vertex_const_iterator;
  using Halfedge_const_iterator =
    typename Arrangement_2::Halfedge_const_iterator;

  using Halfedge_around_vertex_const_circulator =
    typename Arrangement_2::Halfedge_around_vertex_const_circulator;
  using Ccb_halfedge_const_circulator =
    typename Arrangement_2::Ccb_halfedge_const_circulator;
  using Outer_ccb_const_iterator =
    typename Arrangement_2::Outer_ccb_const_iterator;
  using Inner_ccb_const_iterator =
    typename Arrangement_2::Inner_ccb_const_iterator;
  using Isolated_vertex_const_iterator =
    typename Arrangement_2::Isolated_vertex_const_iterator;

  using Point_2 = typename Arrangement_2::Point_2;
  using X_monotone_curve_2 = typename Arrangement_2::X_monotone_curve_2;

  using Result = Arr_point_location_result<Arrangement_2>;
  using Result_type = typename Result::Type;

  // Support cpp11::result_of
  using result_type = Result_type;

private:
  using Gt2 = Geometry_traits_2;
  using Left_side_category =
    typename internal::Arr_complete_left_side_category<Gt2>::Category;
  using Right_side_category =
    typename internal::Arr_complete_right_side_category<Gt2>::Category;
  using Left_or_right_sides_category =
    typename Arr_two_sides_category<Left_side_category,
                                    Right_side_category>::result;

protected:
  using Traits_adaptor_2 = Arr_traits_basic_adaptor_2<Geometry_traits_2>;

  /*! \struct Less_halfedge_handle
   * Used to sort handles.
   */
  struct Less_halfedge_handle {
    bool operator()(Halfedge_const_handle h1, Halfedge_const_handle h2) const
    { return (&(*h1) < &(*h2)); }
  };

  typedef std::set<Halfedge_const_handle, Less_halfedge_handle> Halfedge_set;

  // Data members:
  const Arrangement_2* p_arr;        // The associated arrangement.
  const Traits_adaptor_2* m_traits;  // Its associated traits object.
  Generator* lm_gen;                 // The associated landmark generator.
  bool own_gen;                      // Indicates whether the generator
                                     // has been locally allocated.

  template<typename T>
  Result_type make_result(T t) const { return Result::make_result(t); }
  inline Result_type default_result() const { return Result::default_result(); }

public:
  /*! constructs default. */
  Arr_landmarks_point_location() :
    p_arr(nullptr),
    m_traits(nullptr),
    lm_gen(nullptr),
    own_gen(false)
  {}

  /*! constructs given an arrangement only. */
  Arr_landmarks_point_location(const Arrangement_2& arr) :
    p_arr(&arr),
    m_traits(static_cast<const Traits_adaptor_2*>(p_arr->geometry_traits())),
    lm_gen(new Generator(arr)),         // allocate the landmarks generator.
    own_gen(true)
  {}

  /*! constructs given an arrangement, and landmarks generator. */
  Arr_landmarks_point_location(const Arrangement_2& arr, Generator* gen) :
    p_arr(&arr),
    m_traits(static_cast<const Traits_adaptor_2*>(p_arr->geometry_traits())),
    lm_gen(gen),
    own_gen(false)
  {}

  /*! destructs. */
  ~Arr_landmarks_point_location() {
    if (own_gen) {
      delete lm_gen;
      lm_gen = nullptr;
    }
  }

 /*! attaches an arrangement object (and a generator, if supplied). */
  void attach(const Arrangement_2& arr, Generator* gen = nullptr) {
    // Keep a pointer to the associated arrangement.
    p_arr = &arr;
    m_traits = static_cast<const Traits_adaptor_2*>(p_arr->geometry_traits());

    // Update the landmarks generator.
    if (gen != nullptr) {
      // In case a generator is given, keep a pointer to it.
      CGAL_assertion(lm_gen == nullptr);
      lm_gen = gen;
      own_gen = false;
    }
    else if (lm_gen != nullptr) {
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

  /*! detaches the instance from the arrangement object. */
  void detach() {
    p_arr = nullptr;
    m_traits = nullptr;

    CGAL_assertion(lm_gen != nullptr);
    if (lm_gen)
      lm_gen->detach();
  }

  /*! locates the arrangement feature containing the given point.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type locate(const Point_2& p) const;

protected:
  /*! walks from the given vertex to the query point.
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

  /*! locates an edge around a given vertex that is the predecessor of the
   * curve connecting the vertex to the query point in a clockwise order.
   * \param vh The vertex.
   * \param p The query point.
   * \param new_vertex Output: Whether a closer vertex to p was found.
   * \return The desired object (a halfedge handle or a vertex handle).
   */
  result_type _find_face_around_vertex(Vertex_const_handle vh,
                                       const Point_2& p,
                                       bool& new_vertex) const;

  /*! walks from a point on a given halfedge to the query point.
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
  /*! handles the arrangement curve contained in the segment
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

  /*! walks from a point in a face to the query point.
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

  /*! finds a halfedge on the given CCB that intersects the given x-monotone
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
   *         function returns an invalid halfedge handle).
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

  /*! returns the halfedge that contains the query point.
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

  /*! checks whether the given curve intersects a simple segment, which connects
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

  //!
  template <typename T>
  std::pair<X_monotone_curve_2, Comparison_result>
  construct_segment(const Point_2& p, const Point_2& q, T const& traits,
                    ...) const {
    X_monotone_curve_2 seg = traits.construct_x_monotone_curve_2_object()(p, q);
    Comparison_result res = traits.compare_xy_2_object()(p, q);
    return std::make_pair(seg, res);
  }

  //*!
  template <typename T, typename = typename T::Compare_endpoints_xy_2>
  std::pair<X_monotone_curve_2, Comparison_result>
  construct_segment(const Point_2& p, const Point_2& q, T const& traits,
                    int) const {
    X_monotone_curve_2 seg = traits.construct_x_monotone_curve_2_object()(p, q);
    Comparison_result res = traits.compare_endpoints_xy_2_object()(seg);
    return std::make_pair(seg, res);
  }

  /*! Determines whether the $x$-coordinates of two points are equal.
   */
  bool equal_x_2(const Point_2& p, const Point_2& q,
                 Arr_all_sides_oblivious_tag) const
  { return (m_traits->compare_x_2_object()(p, q) == EQUAL); }

  /*! Determines whether the $x$-coordinates of two points are equal.
   */
  bool equal_x_2(const Point_2& p, const Point_2& q,
                 Arr_has_identified_side_tag) const {
    auto is_on_y_identification = m_traits->is_on_y_identification_2_object();
    if (is_on_y_identification(p)) {
      return is_on_y_identification(q);
    }
    if (is_on_y_identification(q)) return false;
    return (m_traits->compare_x_2_object()(p, q) == EQUAL);
  }

  /*! Determines whether the $x$-coordinates of two points are equal.
   */
  bool equal_x_2(const Point_2& p, const Point_2& q,
                 Arr_boundary_cond_tag) const {
    auto param_space_in_x = m_traits->parameter_space_in_x_2_object();
    switch (param_space_in_x(p)) {
     case ARR_LEFT_BOUNDARY: return (param_space_in_x(q) == ARR_LEFT_BOUNDARY);
     case ARR_RIGHT_BOUNDARY: return (param_space_in_x(q) == ARR_LEFT_BOUNDARY);
     case ARR_INTERIOR: return (m_traits->compare_x_2_object()(p, q) == EQUAL);
     default: CGAL_error();
    }
    CGAL_error();
    return false;
  }
};

} // namespace CGAL

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Arr_landmarks_pl_impl.h>

#include <CGAL/enable_warnings.h>

#endif
