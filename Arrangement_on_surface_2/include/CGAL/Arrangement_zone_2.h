// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ron Wein          <wein@post.tau.ac.il>
//            Efi Fogel         <efif@post.tau.ac.il>
//            (based on old version by Eyal Flato)

#ifndef CGAL_ARRANGEMENT_ZONE_2_H
#define CGAL_ARRANGEMENT_ZONE_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>


/*! \file
 * Defintion of the Arrangement_zone_2 class.
 */

#include <boost/mpl/assert.hpp>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_accessor.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Arr_point_location_result.h>

#include <list>
#include <map>
#include <set>

namespace CGAL {

/*! \class
 * A class for computing the zone of a given $x$-monotone curve in a given
 * arrangement.
 * The arrangement parameter corresponds to the underlying arrangement, and
 * the zone-visitor parameter corresponds to a visitor class which is capable
 * of receiving notifications on the arrangment features the query curve
 * traverses. The visitor has to support the following functions:
 * - init(), for initializing the visitor with a given arrangement.
 * - found_subcurve(), called when a non-intersecting x-monotone curve is
 *                     computed and located in the arrangement.
 * - found_overlap(), called when an x-monotone curve overlaps an existing
 *                    halfedge in the arrangement.
 * Both the second and the third functions return pair<Halfedge_handle, bool>,
 * where the halfedge handle corresponds to the halfedge created or modified
 * by the visitor (if valid), and the Boolean value indicates whether we
 * should halt the zone-computation process.
 */
template <typename Arrangement_, typename ZoneVisitor_>
class Arrangement_zone_2 {
public:
  typedef Arrangement_                                   Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2      Geometry_traits_2;
  typedef typename Arrangement_2::Topology_traits        Topology_traits;

protected:
  typedef Arr_traits_adaptor_2<Geometry_traits_2>        Traits_adaptor_2;

  typedef typename Traits_adaptor_2::Left_side_category   Left_side_category;
  typedef typename Traits_adaptor_2::Bottom_side_category Bottom_side_category;
  typedef typename Traits_adaptor_2::Top_side_category    Top_side_category;
  typedef typename Traits_adaptor_2::Right_side_category  Right_side_category;

  BOOST_MPL_ASSERT
  ((typename Arr_sane_identified_tagging<Left_side_category,
                                         Bottom_side_category,
                                         Top_side_category,
                                         Right_side_category>::result));

public:
  typedef ZoneVisitor_                                   Visitor;

  typedef typename Arrangement_2::Vertex_handle          Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle        Halfedge_handle;
  typedef typename Arrangement_2::Face_handle            Face_handle;

  typedef std::pair<Halfedge_handle, bool>               Visitor_result;

  typedef typename Geometry_traits_2::Point_2            Point_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Geometry_traits_2::Multiplicity       Multiplicity;

protected:
  typedef typename Arr_are_all_sides_oblivious_tag<Left_side_category,
                                                   Bottom_side_category,
                                                   Top_side_category,
                                                   Right_side_category>::result
    Are_all_sides_oblivious_category;

  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;

  typedef typename Arrangement_2::Ccb_halfedge_circulator
    Ccb_halfedge_circulator;

  // Types used for caching intersection points:
  typedef std::pair<Point_2, Multiplicity>              Intersection_point;
  typedef boost::variant<Intersection_point, X_monotone_curve_2>
                                                        Intersection_result;
  typedef boost::optional<Intersection_result>          Optional_intersection;
  typedef std::list<Intersection_result>                Intersect_list;
  typedef std::map<const X_monotone_curve_2*, Intersect_list>
                                                        Intersect_map;
  typedef typename Intersect_map::iterator              Intersect_map_iterator;

  typedef std::set<const X_monotone_curve_2*>           Curves_set;
  typedef typename Curves_set::iterator                 Curves_set_iterator;

  typedef Arr_point_location_result<Arrangement_2>      Pl_result;
  typedef typename Pl_result::Type                      Pl_result_type;

  // Data members:
  Arrangement_2& m_arr;                 // The associated arrangement.
  const Traits_adaptor_2* m_geom_traits; // Its associated geometry traits.
  Arr_accessor<Arrangement_2> m_arr_access; // An accessor for the arrangement.

  Visitor* m_visitor;                   // The zone visitor.

  Intersect_map m_inter_map;            // Stores all computed intersections.

  const Vertex_handle m_invalid_v;      // An invalid vertex handle.
  const Halfedge_handle m_invalid_he;   // An invalid halfedge handle.

  X_monotone_curve_2 m_cv;              // The current portion of the
                                        // inserted curve.
  Pl_result_type m_obj;                 // The location of the left endpoint.
  bool m_has_left_pt;                   // Is the left end of the curve bounded.
  bool m_left_on_boundary;              // Is the left point on the boundary.
  Point_2 m_left_pt;                    // Its current left endpoint.
  bool m_has_right_pt;                  // Is the right end of the curve bounded.
  bool m_right_on_boundary;             // Is the right point on the boundary.
  Point_2 m_right_pt;                   // Its right endpoint (if bounded).

  Vertex_handle m_left_v;               // The arrangement vertex associated
                                        // with the current left_pt (if any).
  Halfedge_handle m_left_he;            // If left_v is valid, left_he is the
                                        // predecessor for cv around this
                                        // vertex. Otherwise, if it is valid,
                                        // it is the halfedge that contains
                                        // the left endpoint it its interior.

  Vertex_handle m_right_v;              // The arrangement vertex associated
                                        // with the current m_right_pt (if any).
  Halfedge_handle m_right_he;           // If m_right_v is valid, left_he is the
                                        // predecessor for cv around this
                                        // vertex. Otherwise, if it is valid,
                                        // it is the halfedge that contains
                                        // the right endpoint it its interior.

  Point_2 m_intersect_p;                // The next intersection point.
  unsigned int m_ip_multiplicity;       // Its multiplicity
                                        // (0 in case of an overlap).
  bool m_found_intersect;               // An intersection has been found.
                                        // (or an overlap).
  X_monotone_curve_2 m_overlap_cv;      // The currently discovered overlap.
  bool m_found_overlap;                 // An overlap has been found.
  bool m_found_iso_vert;                // Check whether an isolated vertex
                                        // induces the next intersection.
  Vertex_handle m_intersect_v;          // The vertex that intersects cv.
  Halfedge_handle m_intersect_he;       // The halfedge that intersects cv
                                        // (or overlaps it).

  X_monotone_curve_2 m_sub_cv1;         // Auxiliary variable (for curve split).
  X_monotone_curve_2 m_sub_cv2;         // Auxiliary variable (for curve split).

public:
  /*! Constructor.
   * \param _arr The arrangement for which we compute the zone.
   * \param _visitor A pointer to a zone-visitor object.
   */
  Arrangement_zone_2(Arrangement_2& arr, Visitor* visitor) :
    m_arr(arr),
    m_arr_access(arr),
    m_visitor(visitor),
    m_invalid_v(),
    m_invalid_he()
  {
    m_geom_traits = static_cast<const Traits_adaptor_2*>(arr.geometry_traits());
    CGAL_assertion(visitor != nullptr);

    // Initialize the visitor.
    visitor->init(&arr);
  }

  /*! Initialize the zone-computation process with a given curve.
   * \param _cv The query curve.
   * \param pl A point-location object associated with the arrangement.
   */
  template <typename PointLocation>
  void init(const X_monotone_curve_2& cv, const PointLocation& pl)
  {
    // Set the curve and check whether its left end has boundary conditions.
    m_cv = cv;

    const Arr_parameter_space bx1 =
      m_geom_traits->parameter_space_in_x_2_object()(m_cv, ARR_MIN_END);
    const Arr_parameter_space by1 =
      m_geom_traits->parameter_space_in_y_2_object()(m_cv, ARR_MIN_END);

    if (bx1 == ARR_INTERIOR && by1 == ARR_INTERIOR) {
      // The curve has a finite left endpoint with no boundary conditions:
      // locate it in the arrangement.
      m_has_left_pt = true;
      m_left_on_boundary = (bx1 != ARR_INTERIOR || by1 != ARR_INTERIOR);
      m_left_pt = m_geom_traits->construct_min_vertex_2_object()(m_cv);

      m_obj = pl.locate(m_left_pt);
    }
    else {
      // The left end of the curve has boundary conditions: use the topology
      // traits use the arrangement accessor to locate it.
      // Note that if the curve-end is unbounded, m_left_pt does not exist.
      // Note that if the curve-end is unbounded, m_left_pt does not exist.
      m_has_left_pt = m_geom_traits->is_closed_2_object()(m_cv, ARR_MIN_END);
      m_left_on_boundary = true;
      if (m_has_left_pt)
        m_left_pt = m_geom_traits->construct_min_vertex_2_object()(m_cv);
      m_obj = m_arr_access.locate_curve_end(m_cv, ARR_MIN_END, bx1, by1);
    }

    // Check the boundary conditions of th right curve end.
    if (m_geom_traits->is_closed_2_object()(m_cv, ARR_MAX_END)) {
      const Arr_parameter_space bx2 =
        m_geom_traits->parameter_space_in_x_2_object()(m_cv, ARR_MAX_END);
      const Arr_parameter_space by2 =
        m_geom_traits->parameter_space_in_y_2_object()(m_cv, ARR_MAX_END);

      // The right endpoint is valid.
      m_has_right_pt = true;
      m_right_pt = m_geom_traits->construct_max_vertex_2_object()(m_cv);
      m_right_on_boundary = (bx2 != ARR_INTERIOR) || (by2 != ARR_INTERIOR);
    }
    else {
      // The right end of the curve lies at infinity.
      m_has_right_pt = false;
      m_right_on_boundary = true;
    }
  }

  /*! Initialize the zone-computation process with a given curve and an object
   * that wraps the location of the curve's left end.
   * \param cv The query curve.
   * \param obj An object that represents the location of the left end of the
   *            curve.
   */
  void init_with_hint(const X_monotone_curve_2& cv, Pl_result_type obj);

  /*! Compute the zone of the given curve and issue the apporpriate
   * notifications for the visitor.
   */
  void compute_zone();

private:
  /*! Check whether two curves with a common endpoint overlap.
   * \pre p == min_point(cv1)
   * \pre p == min_point(cv2)
   * \todo move this function to a more accessible place so that it can be reused
   */
  bool do_overlap(const X_monotone_curve_2& cv1, const X_monotone_curve_2& cv2,
                  const Point_2& p) const
  { return do_overlap_impl(cv1, cv2, p, Are_all_sides_oblivious_category()); }

  /*! Check whether two curves with a common min endpoint overlap.
   */
  bool do_overlap_impl(const X_monotone_curve_2& cv1,
                       const X_monotone_curve_2& cv2,
                       const Point_2& p, Arr_all_sides_oblivious_tag) const
  {
    return m_geom_traits->compare_y_at_x_right_2_object()(cv1, cv2, p) == EQUAL;
  }

  /*! Check whether two curves with a common min endpoint overlap.
   */
  bool do_overlap_impl(const X_monotone_curve_2& cv1,
                       const X_monotone_curve_2& cv2,
                       const Point_2& p, Arr_not_all_sides_oblivious_tag) const;

  /* Check whether the given query curve is encountered when rotating the
   * first curve in a clockwise direction around a given point until reaching
   * the second curve.
   * \pre p == min_point(xcv)
   * \pre p == min_point(xcv1)
   * \pre p == min_point(cxv2)
   * \pre xcv_to_right == TRUE
   * \todo move this function to a more accessible place so that it can be reused
   */
  bool is_between_cw(const X_monotone_curve_2& xcv, bool xcv_to_right,
                     const X_monotone_curve_2& xcv1, bool xcv1_to_right,
                     const X_monotone_curve_2& xcv2, bool xcv2_to_right,
                     const Point_2& p,
                     bool& xcv_equal_xcv1, bool& xcv_equal_xcv2) const
  {
    return is_between_cw_impl(xcv, xcv_to_right,
                              xcv1, xcv1_to_right,
                              xcv2, xcv2_to_right,
                              p, xcv_equal_xcv1, xcv_equal_xcv2,
                              Are_all_sides_oblivious_category());
  }

  /* Check whether the given query curve is encountered when rotating the
   * first curve in a clockwise direction around a given point until reaching
   * the second curve.
   */
  bool is_between_cw_impl(const X_monotone_curve_2& xcv, bool xcv_to_right,
                          const X_monotone_curve_2& xcv1, bool xcv1_to_right,
                          const X_monotone_curve_2& xcv2, bool xcv2_to_right,
                          const Point_2& p,
                          bool& xcv_equal_xcv1, bool& xcv_equal_xcv2,
                          Arr_all_sides_oblivious_tag) const
  {
    return m_geom_traits->is_between_cw_2_object()(xcv, xcv_to_right,
                                                   xcv1, xcv1_to_right,
                                                   xcv2, xcv2_to_right,
                                                   p,
                                                   xcv_equal_xcv1,
                                                   xcv_equal_xcv2);
  }

  /* Check whether the given query curve is encountered when rotating the
   * first curve in a clockwise direction around a given point until reaching
   * the second curve.
   */
  bool is_between_cw_impl(const X_monotone_curve_2& xcv, bool xcv_to_right,
                          const X_monotone_curve_2& xcv1, bool xcv1_to_right,
                          const X_monotone_curve_2& xcv2, bool xcv2_to_right,
                          const Point_2& p,
                          bool& xcv_equal_xcv1, bool& xcv_equal_xcv2,
                          Arr_not_all_sides_oblivious_tag) const;

  /*! Find a face containing the query curve m_cv around the given vertex.
   * In case an overlap occurs, sets m_intersect_he to be the overlapping edge.
   * \param v The query vertex.
   * \param he Output: The predecessor of m_cv around the vertex.
   * \return (true) if m_cv overlaps with the curve associated with he;
   *         (false) if there is no overlap.
   */
  bool _find_prev_around_vertex(Vertex_handle v, Halfedge_handle& he);

  /*! Direct the halfedge for the location of the given subcurve around a split
   * point that occurs in the interior of a given edge, when the subcurve lies
   * to the right of the split point.
   * In case of overlaps, it sets also m_found_overlap and m_intersect_he.
   * \param cv_ins The curve to be inserted, whose left endpoint coincides
   *               with the edge to be split.
   * \param cv_left_pt The left endpoint of cv_ins.
   * \param query_he The edge that intersects cv_ins.
   * \pre The left endpoint of cv_ins lies in the interior of the curve
   *      associated with query_he.
   * \return The halfedge whose incident face contains cv_ins
   *         (either query_he or its twin).
   */
  Halfedge_handle
  _direct_intersecting_edge_to_right(const X_monotone_curve_2& cv_ins,
                                     const Point_2& cv_left_pt,
                                     Halfedge_handle query_he);

  /*! Direct the halfedge for the location of the given subcurve around a split
   * point that occurs in the interior of a given edge, when the subcurve lies
   * to the left of the split point.
   * \param cv_ins The curve to be inserted, whose right endpoint coincides
   *               with the edge to be split.
   * \param query_he The edge that intersects cv_ins.
   * \pre The right endpoint of cv_ins lies in the interior of the curve
   *      associated with query_he.
   * \return The halfedge whose incident face contains cv_ins
   *         (either query_he or its twin).
   */
  Halfedge_handle
  _direct_intersecting_edge_to_left(const X_monotone_curve_2& cv_ins,
                                    Halfedge_handle query_he);

  /*! Get the next intersection of m_cv with the given halfedge.
   * \param he A handle to the halfedge.
   * \param skip_first_point Should we skip the first intersection point.
   * \param intersect_on_right_boundary Output: If an intersetion point is
   *                                            computed, marks whether this
   *                                            point coincides with the right
   *                                            curve-end, which lies on the
   *                                            surface boundary.
   * \return An object representing the next intersection: Intersection_point
   *         in case of a simple intersection point, X_monotone_curve_2 in
   *         case of an overlap, and an empty object if there is no
   *         intersection.
   */
  Optional_intersection
  _compute_next_intersection(Halfedge_handle he,
                             bool skip_first_point,
                             bool& intersect_on_right_boundary);

  /*! Remove the next intersection of m_cv with the given halfedge from the map.
   * \param he A handle to the halfedge.
   * \pre The list of intersections with the curve of he has already been
   *      computed, and it is not empty.
   */
  void _remove_next_intersection (Halfedge_handle he);

  /*! Check whether the given point lies completely to the left of the given
   * egde.
   * \param p The point.
   * \param he The halfedge.
   * \pre he is not a fictitious edge.
   * \return Whether p lies entirely to the left of the edge.
   */
  bool _is_to_left(const Point_2& p, Halfedge_handle he) const
  {
    return (_is_to_left_impl(p, he, Are_all_sides_oblivious_category()));
  }

  bool _is_to_left_impl(const Point_2& p, Halfedge_handle he,
                        Arr_all_sides_oblivious_tag) const
  {
    return (((he->direction() == ARR_LEFT_TO_RIGHT) &&
             m_geom_traits->compare_xy_2_object()
             (p, he->source()->point()) == SMALLER) ||
            (he->direction() == ARR_RIGHT_TO_LEFT &&
             m_geom_traits->compare_xy_2_object()
             (p, he->target()->point()) == SMALLER));
  }

  bool _is_to_left_impl(const Point_2& p, Halfedge_handle he,
                        Arr_not_all_sides_oblivious_tag) const;

  /*! Check whether the given point lies completely to the right of the given
   * egde.
   * \param p The point.
   * \param he The halfedge.
   * \pre he is not a fictitious edge.
   * \return Whether p lies entirely to the right of the edge.
   */
  bool _is_to_right(const Point_2& p, Halfedge_handle he) const
  {
    return (_is_to_right_impl(p, he, Are_all_sides_oblivious_category()));
  }

  bool _is_to_right_impl(const Point_2& p, Halfedge_handle he,
                         Arr_all_sides_oblivious_tag) const
  {
    return (((he->direction() == ARR_LEFT_TO_RIGHT) &&
             m_geom_traits->compare_xy_2_object()(p, he->target()->point()) ==
             LARGER) ||
            ((he->direction() == ARR_RIGHT_TO_LEFT) &&
             m_geom_traits->compare_xy_2_object()(p, he->source()->point()) ==
             LARGER));
  }

  bool _is_to_right_impl(const Point_2& p, Halfedge_handle he,
                         Arr_not_all_sides_oblivious_tag) const;

  /*! Compute the (lexicographically) leftmost intersection of the query
   * curve with a given halfedge on the boundary of a face in the arrangement.
   */
  void
  _leftmost_intersection(Ccb_halfedge_circulator he_curr, bool on_boundary,
                         bool& leftmost_on_right_boundary);

  /*! Compute the (lexicographically) leftmost intersection of the query
   * curve with the boundary of a given face in the arrangement.
   * The function computes sets m_intersect_p, m_intersect_he (or alternatively
   * m_overlap_cv and m_intersect_he) and set the flags m_found_intersect and
   * m_found_overlap accordingly.
   * \param face A handle to the face.
   * \param on_boundary Specifies whether the left endpoint of the curve lies
   *                    on the face boundary.
   */
  void _leftmost_intersection_with_face_boundary(Face_handle face,
                                                 bool on_boundary);

  /*! Compute the zone of an x-monotone curve in a given arrangement face.
   * The left endpoint of the curve either lies in the face interior or on
   * the boundary of the face.
   * This function updates m_cv and its left endpoint and also sets m_left_v
   * and m_left_he for the remaining portion of the curve.
   * In case of overlaps, it sets also m_overlap_cv and m_intersect_he.
   * \param face The given face.
   * \param on_boundary Specifies whether the left endpoint of the curve lies
   *                    on the face boundary.
   * \pre If on_boundary is (true) then m_left_he must be valid; if it is
   *      (false), then both m_left_v anf m_left_he must be invalid.
   * \return (true) if we are done with the zone-computation process;
   *         (false) if we still have a remaining portion of m_cv to continue
   *         with.
   */
  bool _zone_in_face(Face_handle face, bool on_boundary);

  /*! Compute the zone of an overlapping subcurve m_overlap_cv of m_cv and the
   * curve currently associated with m_intersect_he.
   * This function updates m_cv and its left endpoint and also sets m_left_v
   * and m_left_he for the remaining portion of the curve.
   * \return (true) if we are done with the zone-computation process;
   *         (false) if we still have a remaining portion of m_cv to continue
   *         with.
   */
  bool _zone_in_overlap();
};

} //namespace CGAL

// The function definitions can be found under:
#include <CGAL/Arrangement_2/Arrangement_zone_2_impl.h>

#include <CGAL/enable_warnings.h>

#endif
