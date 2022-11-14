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
// Author(s): Baruch Zukerman <baruchzu@post.tau.ac.il>
//            Ron Wein <wein@post.tau.ac.il>
//            Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_NO_INTERSECTION_INSERTION_SS_VISITOR_H
#define CGAL_ARR_NO_INTERSECTION_INSERTION_SS_VISITOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_no_intersection_insertion_ss_visitor class-template.
 * This class can be further split into two, where one derives from the other,
 * such that the derived class handles the case of inserting non-intersecting
 * curves into a non-empty arrangement, and the base class handles the case of
 * inserting non-intersecting curves into a empty arrangement.
 */

#include <CGAL/Surface_sweep_2/Arr_construction_ss_visitor.h>
#include <CGAL/Default.h>

namespace CGAL {

/*! \class Arr_no_intersection_insertion_ss_visitor
 *
 * A sweep-line visitor for inserting new curves into an existing arrangement
 * embedded on a surface, where these curves are interior-disjoint from all
 * existing arrangement edges and vertices (so no intersections occur).
 */
template <typename Helper_, typename Visitor_ = Default>
class Arr_no_intersection_insertion_ss_visitor :
  public Arr_construction_ss_visitor<
    Helper_,
    typename Default::Get<Visitor_, Arr_no_intersection_insertion_ss_visitor<
                                      Helper_, Visitor_> >::type>
{
public:
  typedef Helper_                                       Helper;

  typedef typename Helper::Geometry_traits_2            Geometry_traits_2;
  typedef typename Helper::Event                        Event;
  typedef typename Helper::Subcurve                     Subcurve;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Arr_no_intersection_insertion_ss_visitor<Helper, Visitor_>
                                                        Self;
  typedef typename Default::Get<Visitor_, Self>::type   Visitor;
  typedef Arr_construction_ss_visitor<Helper, Visitor>  Base;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

protected:
  typedef typename Subcurve::Status_line_iterator       Status_line_iterator;
  typedef typename Base::Event_subcurve_reverse_iterator
    Event_subcurve_reverse_iterator;

  typedef typename Helper::Arrangement_2                Arrangement_2;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;

public:
  /*! Constructor. */
  Arr_no_intersection_insertion_ss_visitor(Arrangement_2* arr) : Base(arr) {}

  /// \name Sweep-line notifications.
  //@{

  /* A notification issued before the sweep process starts. */
  void before_sweep();

  /*! A notification invoked before the sweep-line starts handling the given
   * event.
   */
  void before_handle_event(Event* event);

  /*! A notification invoked when a new subcurve is created. */
  void add_subcurve(const X_monotone_curve_2& cv, Subcurve* sc);

  /*! Update the event information. */
  void update_event() {}

  void update_event(Event* /* e */, const Point_2& /* end_point */,
                    const X_monotone_curve_2& /* cv */,
                    Arr_curve_end /* cv_end */, bool /* is_new */)
  {}

  void update_event(Event* /* e */, const X_monotone_curve_2& /* cv */,
                    Arr_curve_end /* cv_end */, bool /* is_new */)
  {}

  void update_event(Event* /* e */, Subcurve* /* sc1 */, Subcurve* /* sc2 */,
                    bool /* is_new */)
  {}

  void update_event(Event* /* e */, Subcurve* /* sc1 */) {}

  void update_event(Event* e, const Point_2& pt, bool /* is_new */)
  {
    Vertex_handle invalid_v;
    if (e->point().vertex_handle() == invalid_v)
      e->point().set_vertex_handle(pt.vertex_handle());
  }
  //@}

  /*! Insert the given subcurve in the interior of a face.
   * \param cv The geometric subcurve.
   * \param sc The sweep-line subcurve information.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle
  insert_in_face_interior(const X_monotone_curve_2& cv, Subcurve* sc);

  /*! Insert the given subcurve given its left end-vertex.
   * \param cv The geometric entity.
   * \param prev The predecessor halfedge around the left vertex.
   * \param sc The sweep-line subcurve information.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle
  insert_from_left_vertex(const X_monotone_curve_2& cv, Halfedge_handle he,
                          Subcurve* sc);

  /*! Insert the given subcurve given its right end-vertex.
   * \param cv The geometric entity.
   * \param prev The predecessor halfedge around the right vertex.
   * \param sc The sweep-line subcurve information.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle
  insert_from_right_vertex(const X_monotone_curve_2& cv, Halfedge_handle prev,
                           Subcurve* sc);

  /*! Insert the given subcurve given its two end-vertices.
   * \param cv The geometric subcurve.
   * \param prev1 The predecessor halfedge around the left vertex.
   * \param prev2 The predecessor halfedge around the right vertex.
   * \param sc The sweep-line subcurve information.
   * \param new_face_created Output: Whether a new face has been created.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle insert_at_vertices(const X_monotone_curve_2& cv,
                                             Halfedge_handle prev1,
                                             Halfedge_handle prev2,
                                             Subcurve* sc,
                                             bool& new_face_created);

  /*! Insert an isolated vertex into the arrangement.
   * \param pt The point associated with the vertex.
   * \param iter The location of the corresponding event in the status line.
   * \return A handle to the inserted vertex.
   */
  virtual Vertex_handle insert_isolated_vertex(const Point_2& pt,
                                               Status_line_iterator iter);

  /// \name Edge-split functions (to be overridden by the child visitor).
  //@{

  /*! Check if the halfedge associated with the given subcurve will be split
   * at the given event.
   * In this case there are no splits.
   */
  virtual bool is_split_event(Subcurve* /*sc*/, Event* /*event*/)
  { return false; }

  /*! Split an edge (does nothing here, as there are no intersections).
   */
  virtual Halfedge_handle split_edge(Halfedge_handle /*he*/,
                                     Subcurve* /*sc*/,
                                     const Point_2& /*pt*/)
  { return Halfedge_handle(); }

  /*! Split an edge (does nothing here, as there are no intersections).
   */
  virtual Halfedge_handle split_edge(Halfedge_handle /*he*/,
                                     Subcurve* /*sc*/,
                                     Vertex_handle /*v*/)
  { return Halfedge_handle(); }
  //@}

protected:
  /// \name Auxiliary functions.
  //@{

  /*! Perform the actual insertion.*/
  Halfedge_handle _insert_in_face_interior(const X_monotone_curve_2& cv,
                                           Subcurve* sc);

  /*! Perform the actual insertion.*/
  Halfedge_handle _insert_from_left_vertex(const X_monotone_curve_2& cv,
                                           Halfedge_handle he,
                                           Subcurve* sc);

  /*! Perform the actual insertion.*/
  Halfedge_handle _insert_from_right_vertex(const X_monotone_curve_2& cv,
                                            Halfedge_handle he,
                                            Subcurve* sc);

  /*! Perform the actual insertion.*/
  Halfedge_handle _insert_at_vertices(const X_monotone_curve_2& cv,
                                      Halfedge_handle hhandle,
                                      Halfedge_handle prev,
                                      Subcurve* sc,
                                      bool& new_face_created);

  /*! Locate the face containing the current object in its interior. */
  Face_handle _ray_shoot_up(Status_line_iterator iter);

  /*! Add a new curve. */
  bool add_subcurve_(const X_monotone_curve_2& cv, Subcurve* sc);
  //@}
};

//-----------------------------------------------------------------------------
// Member-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
// Notifies the helper that the sweep process now starts.
template <typename Hlpr, typename Vis>
void Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::before_sweep()
{ this->m_helper.before_sweep(); }

//-----------------------------------------------------------------------------
// A notification invoked before the sweep-line starts handling the given
// event.
//
template <typename Hlpr, typename Vis>
void Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
before_handle_event(Event* event)
{
  // First we notify the helper class on the event.
  this->m_helper.before_handle_event(event);

  const Halfedge_handle invalid_he;
  event->init_subcurve_in_arrangement_flags(event->number_of_right_curves());
  if (! event->has_right_curves()) {
    // Update the event with the highest left halfedge.
    for (auto it = event->left_curves_rbegin(); it != event->left_curves_rend();
         ++it)
    {
      Halfedge_handle he = (*it)->last_curve().halfedge_handle();
      if (he != invalid_he) {
        event->set_halfedge_handle(he->twin());
        return;
      }
    }
  }

  if (! event->has_left_curves()) {
    // Indicates if there's halfedge to the right of the event.
    int i = 0;
    for (auto it = event->right_curves_rbegin();
         it != event->right_curves_rend(); ++it, ++i)
    {
      // Update the event with the highest right halfedge.
      Halfedge_handle he = (*it)->last_curve().halfedge_handle();
      if (he != invalid_he) {
        event->set_subcurve_in_arrangement(i, true);
        if (event->halfedge_handle() == invalid_he)
          event->set_halfedge_handle(he);
      }
    }
    return;
  }

  // The event has left and right curves.
  bool exist_right_halfedge = false;
  int i = 0;
  for (auto it = event->right_curves_rbegin(); it != event->right_curves_rend();
       ++it, ++i)
  {
    auto he = (*it)->last_curve().halfedge_handle();
    if (he == invalid_he) continue;

    exist_right_halfedge = true;
    event->set_subcurve_in_arrangement(i, true);
    if (! is_split_event(*it, event)) {
      // halfedge will not be split.
      event->set_halfedge_handle(he);
    }
    else {
      Vertex_handle invalid_v;
      he = (event->vertex_handle() != invalid_v) ?
        this->split_edge(he, *it, event->vertex_handle()) :
        this->split_edge(he, *it, event->point());

      // 'he' has the same source as the split halfedge.
      event->set_halfedge_handle(he);
      auto& last_curve = const_cast<X_monotone_curve_2&>((*it)->last_curve());
      last_curve.set_halfedge_handle(he);

      // there cannot be another existing halfedge that need to be split
      // because they are disjoint
      return;
    }
  }

  if (exist_right_halfedge) return;

  // if we have reached here, there are no halfedges to the right of
  // the event, but still can be on the left of the event
  for (auto it = event->left_curves_rbegin();
       it != event->left_curves_rend(); ++it)
  {
    auto he = (*it)->last_curve().halfedge_handle();
    if (he != invalid_he) {
      event->set_halfedge_handle(he->twin());
      return;
    }
  }
}

//-----------------------------------------------------------------------------
// Add a new curve.
//
template <typename Hlpr, typename Vis>
bool Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
add_subcurve_(const X_monotone_curve_2& cv, Subcurve* sc)
{
  const Halfedge_handle invalid_he;
  if (cv.halfedge_handle() != invalid_he) return false;
  // Insert the curve into the arrangement
  Base::add_subcurve(cv, sc);
  return true;
}

//-----------------------------------------------------------------------------
// A notification invoked when a new subcurve is created.
//
template <typename Hlpr, typename Vis>
void Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
add_subcurve(const X_monotone_curve_2& cv, Subcurve* sc)
{
  if (add_subcurve_(cv, sc)) return;

  Halfedge_handle next_ccw_he =
    this->current_event()->halfedge_handle()->next()->twin();
  this->current_event()->set_halfedge_handle(next_ccw_he);
}

//-----------------------------------------------------------------------------
// Insert the given subcurve in the interior of an arrangement face.
//
template <typename Hlpr, typename Vis>
typename Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
insert_in_face_interior(const X_monotone_curve_2& cv, Subcurve* sc)
{
  Event* last_event = this->last_event_on_subcurve(sc);
  Vertex_handle last_v = last_event->point().vertex_handle();
  Vertex_handle curr_v = this->current_event()->point().vertex_handle();
  Vertex_handle null_v;
  if (last_v == null_v && curr_v == null_v)
    return (this->_insert_in_face_interior(cv, sc));
  if (last_v == null_v && curr_v != null_v) {
    Halfedge_handle he =
      this->m_arr->insert_from_right_vertex(cv.base(), curr_v);
    return he->twin();
  }
  if (last_v != null_v && curr_v == null_v)
    return (this->m_arr->insert_from_left_vertex(cv.base(), last_v));
  CGAL_assertion(last_v != null_v && curr_v != null_v);
  return (this->m_arr->insert_at_vertices(cv.base(), last_v, curr_v));
}

//-----------------------------------------------------------------------------
// Insert the given subcurve from a vertex that corresponds to its left end.
//
template <typename Hlpr, typename Vis>
typename Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
insert_from_left_vertex(const X_monotone_curve_2& cv, Halfedge_handle he,
                        Subcurve* sc)
{
  Vertex_handle curr_v = this->current_event()->point().vertex_handle();
  if (curr_v != Vertex_handle())
    return (this->m_arr->insert_at_vertices(cv.base(), he, curr_v));
  return (_insert_from_left_vertex(cv, he, sc));
}

//-----------------------------------------------------------------------------
// Insert the given subcurve from a vertex that corresponds to its right end.
//
template <typename Hlpr, typename Vis>
typename Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
insert_from_right_vertex(const X_monotone_curve_2& cv, Halfedge_handle he,
                         Subcurve* sc)
{
  Event* last_event = this->last_event_on_subcurve(sc);
  Vertex_handle last_v = last_event->point().vertex_handle();
  if (last_v != Vertex_handle())
    return (this->m_arr->insert_at_vertices(cv.base(), he, last_v));
  return (_insert_from_right_vertex(cv, he, sc));
}

//-----------------------------------------------------------------------------
// Insert the given subcurve using its two end-vertices.
//
template <typename Hlpr, typename Vis>
typename Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
insert_at_vertices(const X_monotone_curve_2& cv,
                   Halfedge_handle prev1, Halfedge_handle prev2,
                   Subcurve* sc, bool &new_face_created)
{ return (_insert_at_vertices(cv, prev1, prev2, sc, new_face_created)); }

//-----------------------------------------------------------------------------
// Insert an isolated vertex into the arrangement.
//
template <typename Hlpr, typename Vis>
typename Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::Vertex_handle
Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
insert_isolated_vertex(const Point_2& pt, Status_line_iterator iter)
{
  // If the isolated vertex is already at the arrangement, return:
  if (pt.vertex_handle() != Vertex_handle()) return Vertex_handle();

  // Look up and insert the isolated vertex in the incident face of the
  // halfedge we see.
  Face_handle f = _ray_shoot_up(iter);
  return this->m_arr->insert_in_face_interior(pt.base(), f);
}

//-----------------------------------------------------------------------------
// Perform the actual insertion
//
template <typename Hlpr, typename Vis>
typename Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
_insert_in_face_interior(const X_monotone_curve_2& cv, Subcurve* sc)
{
  // Check if the vertex to be associated with the left end of the curve has
  // already been created.
  Event* last_event = this->last_event_on_subcurve(sc);
  Vertex_handle v1 = last_event->vertex_handle();
  bool create_v1 = false;
  if (v1 == this->m_invalid_vertex) {
    // Mark that we should create the vertex v1 later on (if we created it
    // now, and ended up calling _insert_from_right_vertex(), this vertex
    // would be constructed twice!)
    create_v1 = true;
  }
  else if (v1->degree() > 0) {
    // In this case the left vertex v1 is a boundary vertex which already has
    // some incident halfedges. We look for the predecessor halfedge and
    // and insert the curve from this left vertex.
    Arr_parameter_space bx = last_event->parameter_space_in_x();
    Arr_parameter_space by = last_event->parameter_space_in_y();
    CGAL_assertion((bx != ARR_INTERIOR) || (by != ARR_INTERIOR));
    Halfedge_handle l_prev =
      Halfedge_handle
      (this->m_top_traits->locate_around_boundary_vertex(&(*v1), cv.base(),
                                                         ARR_MIN_END, bx, by));
    return (_insert_from_left_vertex(cv, l_prev, sc));
  }

  // Check if the vertex to be associated with the right end of the curve has
  // already been created.
  Event* curr_event = this->current_event();
  Vertex_handle v2 = curr_event->vertex_handle();

  if (v2 == this->m_invalid_vertex) {
    // Create the vertex to be associated with the right end of the curve.
    v2 = this->m_arr_access.create_vertex(curr_event->point().base());
  }
  else if (v2->degree() > 0) {
    // In this case the right vertex v2 is a boundary vertex which already has
    // some incident halfedges. We look for the predecessor halfedge and
    // and insert the curve from this right vertex.
    Arr_parameter_space bx = curr_event->parameter_space_in_x();
    Arr_parameter_space by = curr_event->parameter_space_in_y();
    CGAL_assertion((bx != ARR_INTERIOR) || (by != ARR_INTERIOR));

    Halfedge_handle r_prev =
      Halfedge_handle
      (this->m_top_traits->locate_around_boundary_vertex(&(*v2), cv.base(),
                                                         ARR_MAX_END, bx, by));
    return (_insert_from_right_vertex(cv, r_prev, sc));
  }

  // If necessary, create the vertex to be associated with the left end
  // of the curve.
  if (create_v1)
    v1 = this->m_arr_access.create_vertex(last_event->point().base());

  // Look up and insert the edge in the interior of the incident face of the
  // halfedge we see.
  Face_handle f = _ray_shoot_up(sc->hint());
  return (this->m_arr_access.insert_in_face_interior_ex(f, cv.base(),
                                                        ARR_LEFT_TO_RIGHT,
                                                        v1, v2));
}

//-----------------------------------------------------------------------------
// Perform the actual insertion
//
template <typename Hlpr, typename Vis>
typename Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
_insert_from_left_vertex(const X_monotone_curve_2& cv,
                         Halfedge_handle prev, Subcurve* sc)
{
  // Check if the vertex to be associated with the right end of the curve has
  // already been created.
  Event* curr_event = this->current_event();
  Vertex_handle v = curr_event->vertex_handle();
  if (v == this->m_invalid_vertex) {
    // Create the vertex to be associated with the right end of the curve.
    v = this->m_arr_access.create_vertex(curr_event->point().base());
  }
  else if (v->degree() > 0) {
    // In this case the left vertex v is a boundary vertex which already has
    // some incident halfedges. We look for the predecessor halfedge and
    // and insert the curve from this right vertex.
    Arr_parameter_space bx = curr_event->parameter_space_in_x();
    Arr_parameter_space by = curr_event->parameter_space_in_y();
    CGAL_assertion(bx != ARR_INTERIOR || by != ARR_INTERIOR);
    Halfedge_handle r_prev =
      Halfedge_handle
      (this->m_top_traits->locate_around_boundary_vertex(&(*v), cv.base(),
                                                         ARR_MAX_END, bx, by));
    bool dummy;
    return (_insert_at_vertices(cv, r_prev, prev, sc, dummy));
  }

  // Perform the insertion using the vertex v.
  return (this->m_arr_access.insert_from_vertex_ex(prev, cv.base(),
                                                   ARR_LEFT_TO_RIGHT, v));
}

//-----------------------------------------------------------------------------
// Perform the actual insertion
//
template <typename Hlpr, typename Vis>
typename Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
_insert_from_right_vertex(const X_monotone_curve_2& cv, Halfedge_handle prev,
                          Subcurve* sc)
{
  // Check if the vertex to be associated with the left end of the curve has
  // already been created.
  Event* last_event = this->last_event_on_subcurve(sc);
  Vertex_handle v = last_event->vertex_handle();
  if (v == this->m_invalid_vertex) {
    // Create the vertex to be associated with the left end of the curve.
    v = this->m_arr_access.create_vertex(last_event->point().base());
  }
  else if (v->degree() > 0) {
    // In this case the left vertex v is a boundary vertex which already has
    // some incident halfedges. We look for the predecessor halfedge and
    // and insert the curve between two existing vertices.
    Arr_parameter_space bx = last_event->parameter_space_in_x();
    Arr_parameter_space by = last_event->parameter_space_in_y();
    CGAL_assertion((bx != ARR_INTERIOR) || (by != ARR_INTERIOR));
    Halfedge_handle l_prev =
      Halfedge_handle
      (this->m_top_traits->locate_around_boundary_vertex(&(*v), cv.base(),
                                                         ARR_MIN_END, bx, by));
    bool dummy;
    return (_insert_at_vertices(cv, prev, l_prev, sc, dummy));
  }

  // Perform the insertion using the vertex v.
  return (this->m_arr_access.insert_from_vertex_ex(prev, cv.base(),
                                                   ARR_RIGHT_TO_LEFT, v));
}

//-----------------------------------------------------------------------------
// Perform the actual insertion
//
template <typename Hlpr, typename Vis>
typename Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
_insert_at_vertices(const X_monotone_curve_2& cv,
                    Halfedge_handle prev1, Halfedge_handle prev2,
                    Subcurve* sc, bool& new_face_created)
{
  // Perform the insertion.
  new_face_created = false;
  bool swapped_predecessors = false;
  Halfedge_handle new_he =
    this->m_arr_access.insert_at_vertices_ex(prev1,
                                             cv.base(), ARR_RIGHT_TO_LEFT,
                                             prev2->next(),
                                             new_face_created,
                                             swapped_predecessors,
                                             true /* this allows to swap
                                                   * prev1/prev2 which is done
                                                   * by checking local minima,
                                                   * which is cheaper than
                                                   * comparing lengths */
                                             );

  // Notify the helper on the creation of the new halfedge.
  // Note that we do this before trying to relocate holes in the new
  // face (if one is created to match the corresponding call from the
  // construction visitor).
  this->m_helper.add_subcurve(new_he, sc);

  // In case a new face has been created (pointed by the new halfedge we
  // obtained), we have to examine the holes and isolated vertices in the
  // existing face (pointed by the twin halfedge) and move the relevant
  // holes and isolated vertices into the new face.
  if (new_face_created) this->m_arr_access.relocate_in_new_face(new_he);

  // Return a handle to the new halfedge directed from prev1's target to
  // prev2's target. Note that this may be the twin halfedge of the one
  // returned by _insert_at_vertices();
  if (swapped_predecessors) new_he = new_he->twin();
  return new_he;
}

//-----------------------------------------------------------------------------
// Locate the face containing the current object in its interior.
//
template <typename Hlpr, typename Vis>
typename Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::Face_handle
Arr_no_intersection_insertion_ss_visitor<Hlpr, Vis>::
_ray_shoot_up(Status_line_iterator iter)
{
  // Go up the status line and try to locate a curve which is associated
  // with a valid arrangement halfedge.
  const Halfedge_handle invalid_he;
  for (; iter != this->status_line_end(); ++iter) {
    Halfedge_handle he = (*iter)->last_curve().halfedge_handle();
    // Return the incident face of the halfedge if found.
    if (he != invalid_he) return (he->face());
  }

  // If we reached here, there is no arrangement halfedge above the given
  // subcurve. We therefore return the top face, as given by the helper class.
  return (this->m_helper.top_face());
}

} // namespace CGAL

#endif
