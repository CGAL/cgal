// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_OVERLAY_SL_VISITOR_H
#define CGAL_ARR_OVERLAY_SL_VISITOR_H

/*! \file
 * Definition of the Arr_overlay_sl_visitor class-template.
 */

#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/unordered_map.hpp>
#include <boost/variant/apply_visitor.hpp>

#include <CGAL/Arr_tags.h>
#include <CGAL/Sweep_line_2/Arr_construction_sl_visitor.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

/*! \class Arr_overlay_sl_visitor
 * A sweep-line visitor for overlaying a "red" arrangement and a "blue"
 * arrangement, creating a result arrangement. All three arrangements are
 * embedded on the same type of surface and use the same geometry traits.
 */
template <typename OverlayHelper_, typename OverlayTraits_>
class Arr_overlay_sl_visitor : public
  Arr_construction_sl_visitor<typename OverlayHelper_::Construction_helper>
{
public:
  typedef OverlayHelper_                                   Overlay_helper;
  typedef OverlayTraits_                                   Overlay_traits;

  typedef typename Overlay_helper::Traits_2                Traits_2;
  typedef typename Overlay_helper::Event                   Event;
  typedef typename Overlay_helper::Subcurve                Subcurve;

  typedef typename Traits_2::X_monotone_curve_2            X_monotone_curve_2;
  typedef typename Traits_2::Point_2                       Point_2;

  // The input arrangements (the "red" and the "blue" one):
  typedef typename Overlay_helper::Arrangement_red_2       Arrangement_red_2;
  typedef typename Arrangement_red_2::Halfedge_const_handle
                                                           Halfedge_handle_red;
  typedef typename Arrangement_red_2::Face_const_handle    Face_handle_red;
  typedef typename Arrangement_red_2::Vertex_const_handle  Vertex_handle_red;

  typedef typename Overlay_helper::Arrangement_blue_2      Arrangement_blue_2;
  typedef typename Arrangement_blue_2::Halfedge_const_handle
                                                           Halfedge_handle_blue;
  typedef typename Arrangement_blue_2::Face_const_handle   Face_handle_blue;
  typedef typename Arrangement_blue_2::Vertex_const_handle Vertex_handle_blue;

  // The resulting arrangement:
  typedef typename Overlay_helper::Arrangement_2           Arrangement_2;
  typedef typename Arrangement_2::Halfedge_handle          Halfedge_handle;
  typedef typename Arrangement_2::Face_handle              Face_handle;
  typedef typename Arrangement_2::Vertex_handle            Vertex_handle;
  typedef typename Arrangement_2::Ccb_halfedge_circulator
    Ccb_halfedge_circulator;
  typedef typename Arrangement_2::Outer_ccb_iterator       Outer_ccb_iterator;

  // The base construction visitor:
  typedef typename Overlay_helper::Construction_helper     Construction_helper;
  typedef Arr_construction_sl_visitor<Construction_helper> Base;

  typedef typename Base::Event_subcurve_iterator
    Event_subcurve_iterator;
  typedef typename Base::Event_subcurve_reverse_iterator
    Event_subcurve_reverse_iterator;
  typedef typename Base::Status_line_iterator             Status_line_iterator;

protected:
  typedef typename Traits_2::Cell_handle_red              Cell_handle_red;
  typedef typename Traits_2::Optional_cell_red            Optional_cell_red;
  typedef typename Traits_2::Cell_handle_blue             Cell_handle_blue;
  typedef typename Traits_2::Optional_cell_blue           Optional_cell_blue;

  typedef std::pair<Halfedge_handle_red, Halfedge_handle_blue>
                                                          Halfedge_info;
  typedef Unique_hash_map<Halfedge_handle,Halfedge_info>  Halfedge_map;

  typedef std::pair<Cell_handle_red, Cell_handle_blue>    Handle_info;
  typedef boost::unordered_map<Vertex_handle, Handle_info, Handle_hash_function>
                                                          Vertex_map;

  // Side categoties:
  typedef typename Traits_2::Left_side_category           Left_side_category;
  typedef typename Traits_2::Bottom_side_category         Bottom_side_category;
  typedef typename Traits_2::Top_side_category            Top_side_category;
  typedef typename Traits_2::Right_side_category          Right_side_category;

  typedef typename Arr_has_identified_sides<Left_side_category,
                                            Bottom_side_category>::result
    Has_identified_sides_category;

  // Data members:
  Overlay_traits* m_overlay_traits;     // The overlay traits object.

  Overlay_helper m_overlay_helper;      // The overlay-helper class
                                        // (note that the base class stores
                                        //  an additional helper class used
                                        //  for constructing the result).

  Halfedge_map m_halfedges_map;         // Mapping of each halfedge in the
                                        // result arrangement to the red
                                        // and blue halfedges that induce it.

  Vertex_map m_vertices_map;            // Mapping of boundary vertices in
                                        // the result arrangement to the red
                                        // and blue halfedges that induce it.
public:
  /*! Constructor */
  Arr_overlay_sl_visitor(const Arrangement_red_2* red_arr,
                         const Arrangement_blue_2* blue_arr,
                         Arrangement_2* res_arr,
                         Overlay_traits* overlay_traits):
    Base(res_arr),
    m_overlay_traits(overlay_traits),
    m_overlay_helper(red_arr, blue_arr),
    m_halfedges_map(Halfedge_info(),
                    // Give an initial size for the hash table
                    red_arr->number_of_halfedges() +
                    blue_arr->number_of_halfedges())
  {}

  /*! Destructor */
  virtual ~Arr_overlay_sl_visitor() {}

  /// \name Sweep-line notifications.
  //@{

  /* A notification issued before the sweep process starts. */
  void before_sweep();

  /*!
   * A notification invoked before the sweep-line starts handling the given
   * event.
   */
  void before_handle_event(Event* event);

  /*!
   * A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  bool after_handle_event(Event* event, Status_line_iterator iter, bool flag);

  /*! Update an event that corresponds to a curve endpoint. */
  void update_event(Event* e,
                    const Point_2& end_point,
                    const X_monotone_curve_2& cv,
                    Arr_curve_end cv_end,
                    bool is_new);

   void update_event(Event* /* e */,
                     const X_monotone_curve_2& /* cv */,
                     Arr_curve_end /* cv_end */,
                     bool /* is_new */)
   {}

  /*! Update an event that corresponds to an intersection between curves. */
  void update_event(Event* /* e */,
                    Subcurve* /* c1 */,
                    Subcurve* /* c2 */,
                    bool CGAL_assertion_code(is_new))
  {
    CGAL_assertion(is_new == true);
  }

  /*! Update an event. */
  void update_event(Event* e, Subcurve* sc);

  /*! Update an event. */
  void update_event(Event* e, const Point_2& p, bool is_new);

  /* A notification issued when the sweep process has ended. */
  void after_sweep();
  //@}

  /// \name Insertion functions.
  //@{

  /*!
   * Insert the given subcurve in the interior of a face.
   * \param cv The geometric subcurve.
   * \param sc The sweep-line subcurve information.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle insert_in_face_interior(const X_monotone_curve_2& cv,
                                                  Subcurve* sc);

  /*!
   * Insert the given subcurve given its right end-vertex.
   * \param cv The geometric entity.
   * \param prev The predecessor halfedge around the right vertex.
   * \param sc The sweep-line subcurve information.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle insert_from_right_vertex(const X_monotone_curve_2& cv,
                                                   Halfedge_handle prev,
                                                   Subcurve* sc);

  /*!
   * Insert the given subcurve given its left end-vertex.
   * \param cv The geometric entity.
   * \param prev The predecessor halfedge around the left vertex.
   * \param sc The sweep-line subcurve information.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle insert_from_left_vertex(const X_monotone_curve_2& cv,
                                                  Halfedge_handle prev,
                                                  Subcurve* sc);

 /*!
   * Insert the given subcurve given its two end-vertices.
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

  /*!
   * Insert an isolated vertex into the arrangement.
   * \param pt The point associated with the vertex.
   * \param iter The location of the corresponding event in the status line.
   * \return A handle to the inserted vertex.
   */
  virtual Vertex_handle insert_isolated_vertex(const Point_2& pt,
                                               Status_line_iterator iter);
  //@}

protected:
  /// \name Auxiliary functions.
  //@{

  /*!
   * Map a newly created halfedge in the result arrangement to its originator
   * red and blue halfedges.
   * \param he The new halfedge.
   * \param red_he The originating red halfedge.
   * \param blue_he The originating blue halfedge.
   * \pre he is directed right to left.
   */
  void _map_halfedge_and_twin(Halfedge_handle he,
                              Halfedge_handle_red red_he,
                              Halfedge_handle_blue blue_he);

  /*!
   * Update the boundary vertices map.
   * This function is used when the parameter space has an identified (or
   * contructed) boundary side. We assume that if the parameter space has a
   * contructed boundary side, it also must have an identified boundary side.
   * \param event The event.
   * \param v The vertex.
   * \param tag The tag used for dispatching.
   */
  void _map_boundary_vertices(Event* event, Vertex_handle v,
                              boost::mpl::bool_<true> /* tag */);

  /*!
   * Update the boundary vertices map.
   * This function is used when the parameter space does not have an identified
   * boundary side, and thus, neither it has a contructed boundary side.
   * \param event The event.
   * \param v The vertex.
   * \param tag The tag used for dispatching.
   */
  void _map_boundary_vertices(Event* event, Vertex_handle v,
                              boost::mpl::bool_<false> /* tag */);

  /*!
   * Update a newly created vertex using the overlay traits.
   * This function is used when the parameter space has an identified (or
   * contructed) boundary side. We assume that if the parameter space has a
   * contructed boundary side, it also must have an identified boundary side.
   * \param event The event associated with the new vertex.
   * \param res_v The new vertex in the overlaid arrangement.
   * \param sc The subcurve incident to the event.
   * \param tag The tag used for dispatching.
   */
  void _create_vertex(Event* event, Vertex_handle res_v, Subcurve* sc,
                      boost::mpl::bool_<true> /* tag */);

  /*!
   * Update a newly created vertex using the overlay traits.
   * This function is used when the parameter space does not have an identified
   * boundary side, and thus, neither it has a contructed boundary side.
   * \param event The event associated with the new vertex.
   * \param res_v The new vertex in the overlaid arrangement.
   * \param sc The subcurve incident to the event.
   * \param tag The tag used for dispatching.
   */
  void _create_vertex(Event* event, Vertex_handle res_v, Subcurve* sc,
                      boost::mpl::bool_<false> /* tag */);

  /*!
   * Update a newly created edge using the overlay traits.
   * \param sc The subcurve associated with the new edge.
   * \param res_he One of the new halfedges in the overlaid arrangement.
   */
  void _create_edge(Subcurve* sc, Halfedge_handle res_he);

  /*! A visitor class to facilitate the call to create_vertex(). */
  class Create_vertex_visitor : public boost::static_visitor<> {
  private:
    Overlay_traits* m_overlay_traits;
    Vertex_handle m_vertex_handle;

  public:
    /*! Constructor */
    Create_vertex_visitor(Overlay_traits* overlay_traits,
                          Vertex_handle vertex_handle) :
      m_overlay_traits(overlay_traits),
      m_vertex_handle(vertex_handle)
    {}

    /*! */
    void operator()(const Vertex_handle_red& red_v,
                    const Vertex_handle_blue& blue_v) const
    { m_overlay_traits->create_vertex(red_v, blue_v, m_vertex_handle); }

    /*! */
    void operator()(const Halfedge_handle_red& red_he,
                    const Vertex_handle_blue& blue_v) const
    { m_overlay_traits->create_vertex(red_he, blue_v, m_vertex_handle); }

    /*! */
    void operator()(const Face_handle_red& red_f,
                    const Vertex_handle_blue& blue_v) const
    { m_overlay_traits->create_vertex(red_f, blue_v, m_vertex_handle); }

    /*! */
    void operator()(const Vertex_handle_red& red_v,
                    const Halfedge_handle_blue& blue_he) const
    { m_overlay_traits->create_vertex(red_v, blue_he, m_vertex_handle); }

    /*! */
    void operator()(const Vertex_handle_red& red_v,
                    const Face_handle_blue& blue_f) const
    { m_overlay_traits->create_vertex(red_v, blue_f, m_vertex_handle); }

    /*! */
    void operator()(const Halfedge_handle_red& red_he,
                    const Halfedge_handle_blue& blue_he) const
    { m_overlay_traits->create_vertex(red_he, blue_he, m_vertex_handle); }

    /* The following functions should never be called */
    void operator()(const Halfedge_handle_red& /* red_he */,
                    const Face_handle_blue& /* blue_f */) const
    { CGAL_error(); }
    void operator()(const Face_handle_red& /* red_f */,
                    const Halfedge_handle_blue& /* blue_he */) const
    { CGAL_error(); }
    void operator()(const Face_handle_red& /* red_f */,
                    const Face_handle_blue& /* blue_f */) const
    { CGAL_error(); }

  };
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
template <typename OvlHlpr, typename OvlTr>
void Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::before_sweep()
{
  // Initialize the necessary fields in the base construction visitor.
  // Note that the construction visitor also informs its helper class that
  // the sweep process is about to start.
  Base::before_sweep();

  // Notify the overlay helper that the sweep process now starts.
  m_overlay_helper.before_sweep();
}

//-----------------------------------------------------------------------------
// A notification invoked before the sweep-line starts handling the given
// event.
//
template <typename OvlHlpr, typename OvlTr>
void Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::before_handle_event(Event* event)
{
  // Let the base construction visitor do the work (and also inform its helper
  // class on the event).
  Base::before_handle_event(event);

  // Notify the overlay helper class on the event.
  m_overlay_helper.before_handle_event(event);
}

//-----------------------------------------------------------------------------
// A notification invoked after the sweep-line finishes handling the given
// event.
//
template <typename OvlHlpr, typename OvlTr>
bool Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
after_handle_event(Event* event, Status_line_iterator iter, bool flag)
{
  // Let the base construction visitor handle the event.
  bool res = Base::after_handle_event(event, iter, flag);

  // In case there are no subcurves in the status line above the given event
  // point, we update the top fictitious halfedges for all subcurves incident
  // to this event.
  Event_subcurve_reverse_iterator rev_iter = event->right_curves_rbegin();
  Subcurve* sc_above = NULL;

  if (iter != this->status_line_end())
    sc_above = (*iter);

  if (sc_above == NULL) {
    if (rev_iter != event->right_curves_rend()) {
      if ((*rev_iter)->color() == Traits_2::BLUE)
        (*rev_iter)->set_red_top_face(m_overlay_helper.red_top_face());
      else if ((*rev_iter)->color() == Traits_2::RED)
        (*rev_iter)->set_blue_top_face(m_overlay_helper.blue_top_face());

      (*rev_iter)->set_subcurve_above(NULL);
      sc_above = *rev_iter;
      ++rev_iter;
    }
    else return res;            // Nothing else to do.
  }

  // For each subcurve, try to locate a subcurve lying above it in the status
  // line that has a different color.
  for (; rev_iter != event->right_curves_rend(); ++rev_iter) {
    Subcurve* curr_sc = *rev_iter;

    if (! curr_sc->has_same_color(sc_above))
      curr_sc->set_subcurve_above(sc_above);
    else {
      if (sc_above->subcurve_above() != NULL)
        curr_sc->set_subcurve_above(sc_above->subcurve_above());
      else
        curr_sc->set_top_face(sc_above);
    }

    sc_above = curr_sc;
  }
  return (res);
}

//-----------------------------------------------------------------------------
// Update an event that corresponds to a curve endpoint.
//
template <typename OvlHlpr, typename OvlTr>
void Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
update_event(Event* e,
             const Point_2& end_point,
             const X_monotone_curve_2& /* cv */,
             Arr_curve_end /* cv_end */,
             bool /* is_new */)
{
  // Nothing to do in case of an event at infinity.
  CGAL_assertion(e->is_closed());

  // Update the red and blue objects associated with the point as necessary.
  Point_2& pt = e->point();
  if (pt.is_red_cell_empty()) pt.set_red_cell(end_point.red_cell());
  else if (pt.is_blue_cell_empty()) pt.set_blue_cell(end_point.blue_cell());
}

//-----------------------------------------------------------------------------
// Update an event.
//
template <typename OvlHlpr, typename OvlTr>
void Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::update_event(Event* e,
                                                          Subcurve* sc)
{
  // Update the red and blue halfedges associated with the point as necessary.
  Point_2& pt = e->point();

  if (pt.is_red_cell_empty()) {
    CGAL_assertion(! pt.is_blue_cell_empty());
    CGAL_assertion(sc->color() == Traits_2::RED);

    Halfedge_handle_red red_he = sc->red_halfedge_handle();
    pt.set_red_cell(boost::make_optional(Cell_handle_red(red_he)));
  }
  else if (pt.is_blue_cell_empty()) {
    Halfedge_handle_blue blue_he = sc->blue_halfedge_handle();
    pt.set_blue_cell(boost::make_optional(Cell_handle_blue(blue_he)));
  }
}

//-----------------------------------------------------------------------------
// Update an event.
//
template <typename OvlHlpr, typename OvlTr>
void Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::update_event(Event* e,
                                                          const Point_2& p,
                                                          bool /* is_new */)
{
  // Update the red and blue objects associated with the point as necessary.
  Point_2& pt = e->point();
  if (pt.is_red_cell_empty()) pt.set_red_cell(p.red_cell());
  else if (pt.is_blue_cell_empty()) pt.set_blue_cell(p.blue_cell());
}

//-----------------------------------------------------------------------------
// A notification issued when the sweep process has ended.
//
template <typename OvlHlpr, typename OvlTr>
void Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::after_sweep()
{
  // Notify boundary vertices:
  typename Vertex_map::iterator it;
  for (it = m_vertices_map.begin(); it != m_vertices_map.end(); ++it) {
    const Handle_info& info = (*it).second;
    const Cell_handle_red& red_handle = info.first;
    const Cell_handle_blue& blue_handle = info.second;
    Vertex_handle v = (*it).first;
    Create_vertex_visitor visitor(m_overlay_traits, v);
    boost::apply_visitor(visitor, red_handle, blue_handle);
  }

  // When the sweep-line process is over, the remaining arrangement face
  // (the current top face of the result arrangement) should be updated such
  // that it is the overlay of the two remaining red and blue unbounded face
  // (the current "red" and "blue" top faces, as given by the helper class).
  m_overlay_traits->create_face(m_overlay_helper.red_top_face(),
                                m_overlay_helper.blue_top_face(),
                                this->m_helper.top_face());
}

//-----------------------------------------------------------------------------
// Insert the given subcurve in the interior of an arrangement face.
//
template <typename OvlHlpr, typename OvlTr>
typename Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::Halfedge_handle
Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
insert_in_face_interior(const X_monotone_curve_2& cv, Subcurve* sc)
{
  // Insert the halfedge using the base construction visitor.
  Halfedge_handle new_he = Base::insert_in_face_interior(cv, sc);
  _map_halfedge_and_twin(new_he,
                         cv.red_halfedge_handle(), cv.blue_halfedge_handle());

  // Update the newly created left and right vertices using the overlay traits.
  if (new_he->direction() == ARR_LEFT_TO_RIGHT) {
    _create_vertex(this->last_event_on_subcurve(sc), new_he->source(), sc,
                   Has_identified_sides_category());
    _create_vertex(this->current_event(), new_he->target(), sc,
                   Has_identified_sides_category());
  }
  else {
    _create_vertex(this->last_event_on_subcurve(sc), new_he->target(), sc,
                   Has_identified_sides_category());
    _create_vertex(this->current_event(), new_he->source(), sc,
                   Has_identified_sides_category());
  }

  // Update the newly created edge using the overlay traits.
  _create_edge(sc, new_he);

  return new_he;
}

//-----------------------------------------------------------------------------
// Insert the given subcurve given its left end-vertex.
//
template <typename OvlHlpr, typename OvlTr>
typename Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::Halfedge_handle
Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
insert_from_left_vertex(const X_monotone_curve_2& cv,
                        Halfedge_handle prev,
                        Subcurve* sc)
{
  _map_boundary_vertices(this->last_event_on_subcurve(sc), prev->target(),
                         Has_identified_sides_category());

  // Insert the halfedge using the base construction visitor. Note that the
  // resulting halfedge is directed from left to right.
  Halfedge_handle new_he = Base::insert_from_left_vertex(cv, prev, sc);
  _map_halfedge_and_twin(new_he,
                         cv.red_halfedge_handle(), cv.blue_halfedge_handle());

  // Update the newly created vertex, which is on the right, using the overlay
  // traits.
  Vertex_handle new_v_right = (new_he->direction() == ARR_LEFT_TO_RIGHT) ?
    new_he->target() : new_he->source();
  _create_vertex(this->current_event(), new_v_right, sc,
                 Has_identified_sides_category());

  // Update the newly created edge using the overlay traits.
  _create_edge(sc, new_he);

  return new_he;
}

//-----------------------------------------------------------------------------
// Insert the given subcurve given its right end-vertex.
//
template <typename OvlHlpr, typename OvlTr>
typename Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::Halfedge_handle
Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
insert_from_right_vertex(const X_monotone_curve_2& cv,
                         Halfedge_handle prev,
                         Subcurve* sc)
{
  _map_boundary_vertices(this->current_event(), prev->target(),
                         Has_identified_sides_category());

  // Insert the halfedge using the base construction visitor.
  Halfedge_handle new_he = Base::insert_from_right_vertex(cv, prev, sc);
  _map_halfedge_and_twin(new_he,
                         cv.red_halfedge_handle(), cv.blue_halfedge_handle());

  // Update the newly created vertex, which is on the left, using the overlay
  // traits.
  Vertex_handle new_v_left = (new_he->direction() == ARR_LEFT_TO_RIGHT) ?
    new_he->source() : new_he->target();
  _create_vertex(this->last_event_on_subcurve(sc), new_v_left, sc,
                 Has_identified_sides_category());

  // Update the newly created edge using the overlay traits.
  _create_edge(sc, new_he);

  return new_he;
}

//-----------------------------------------------------------------------------
// Insert the given subcurve given its two end-vertices.
//
template <typename OvlHlpr, typename OvlTr>
typename Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::Halfedge_handle
Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
insert_at_vertices(const X_monotone_curve_2& cv,
                   Halfedge_handle prev1,
                   Halfedge_handle prev2,
                   Subcurve* sc,
                   bool& new_face_created)
{
  // Insert the halfedge using the base construction visitor. Note that the
  // resulting halfedge is always incident to the new face (if one created).
  Halfedge_handle new_he =
    Base::insert_at_vertices(cv, prev1, prev2, sc, new_face_created);
  _map_halfedge_and_twin(new_he,
                         cv.red_halfedge_handle(), cv.blue_halfedge_handle());

  _map_boundary_vertices(this->last_event_on_subcurve(sc),
                         (new_he->direction() == ARR_LEFT_TO_RIGHT) ?
                         new_he->source() : new_he->target(),
                         Has_identified_sides_category());

  _map_boundary_vertices(this->current_event(),
                         (new_he->direction() == ARR_LEFT_TO_RIGHT) ?
                         new_he->target() : new_he->source(),
                         Has_identified_sides_category());

  // Update the newly created edge using the overlay traits.
  _create_edge(sc, new_he);

  // If a new face was created, we have to updated the newly created face
  // using the overlay traits.
  if (new_face_created) {
    // Obtain the new face, which is incident to new_he.
    Face_handle new_face = new_he->face();

    // Traverse the boundary of the new face, and locate halfedge originated
    // by red or by blue halfedges along its boundary.
    // We stop the traversal earlier if we locate a red and a blue halfedge.
    const Halfedge_handle_red   invalid_red_he;
    const Halfedge_handle_blue  invalid_blue_he;

    Halfedge_handle_red         red_he;
    Halfedge_handle_blue        blue_he;

    CGAL_assertion(new_face->number_of_outer_ccbs() > 0);
    Outer_ccb_iterator occb_it = new_face->outer_ccbs_begin();
    // msvc CL requires the breakdown to the following 2 statements:
    Ccb_halfedge_circulator     ccb_first = *occb_it;
    Ccb_halfedge_circulator     ccb_circ = ccb_first;

    do {
      // Get the current halfedge on the face boundary and obtain its
      // originating halfedges.
      Halfedge_handle he = ccb_circ;

      if (! m_halfedges_map.is_defined(he)) {
        // The mapping is not available for fictitious halfedges ...
        ++ccb_circ;
        continue;
      }

      // The halfedge info is a pair of red and blue halfedges, one of them
      // may be invalid.
      const Halfedge_info& he_info = m_halfedges_map[he];

      if (he_info.first != invalid_red_he) {
        red_he = he_info.first;
        if (blue_he != invalid_blue_he) break;
      }

      if (he_info.second != invalid_blue_he) {
        blue_he = he_info.second;
        if (red_he != invalid_red_he) break;
      }

      ++ccb_circ;
    } while (ccb_circ != ccb_first);

    // Determine the relation between the red and blue face that originate
    // our new overlay face and update it accordingly.
    Face_handle_red   red_face;
    Face_handle_blue  blue_face;

    if ((red_he != invalid_red_he) && (blue_he != invalid_blue_he)) {
      // The red face and the blue face intersect (or overlap).
      red_face = red_he->face();
      blue_face = blue_he->face();
    }
    else if (red_he != invalid_red_he) {
      // The new overlay face is an entirely red face contained inside a blue
      // face. We have to find the identity of this containing blue face.
      Subcurve* sc_above = sc->subcurve_above();
      red_face = red_he->face();
      blue_face = (sc_above != NULL) ?
        sc_above->blue_halfedge_handle()->face() : sc->blue_top_face();
    }
    else {
      CGAL_assertion(blue_he != invalid_blue_he);

      // The new overlay face is an entirely blue face contained inside a red
      // face. We have to find the identity of this containing red face.
      Subcurve* sc_above = sc->subcurve_above();
      blue_face = blue_he->face();
      red_face = (sc_above != NULL) ?
        sc_above->red_halfedge_handle()->face() : sc->red_top_face();
    }

    // Use the overlay traits to update the information of the newly created
    // face.
    m_overlay_traits->create_face(red_face, blue_face, new_face);
  }

  return new_he;
}

//-----------------------------------------------------------------------------
// Insert an isolated vertex into the arrangement.
//
template <typename OvlHlpr, typename OvlTr>
typename Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::Vertex_handle
Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
insert_isolated_vertex(const Point_2& pt,
                       Status_line_iterator iter)
{
  // Insert the isolated vertex using the base construction visitor.
  Vertex_handle new_v = Base::insert_isolated_vertex(pt, iter);

  // Get the red and blue objects that induced the isolated point.
  // Note that as this point is isolated, both objects must be vertices
  // (in case they are not empty).
  const Vertex_handle_red* red_vertex = pt.red_vertex_handle();
  const Vertex_handle_blue* blue_vertex = pt.blue_vertex_handle();

  CGAL_assertion(red_vertex || blue_vertex);

  if (red_vertex && blue_vertex) {
    // The vertex is created by two coincident isolated vertices.
    // Use the overlay traits to update the newly created isolated vertex.
    m_overlay_traits->create_vertex(*red_vertex, *blue_vertex, new_v);
  }
  else if (! red_vertex) {
    // We have an isolated blue vertex inside a red face.
    Face_handle_red red_face;

    // Obtain the red face containing the isolated vertex.
    if (iter == this->status_line_end()) {
      // There is nothing above the vertex - use the current red top face.
      red_face = m_overlay_helper.red_top_face();
    }
    else {
      // Go up the status line until we find a red halfedge.
      // Note that the subcurve above is always of an opposite color. However,
      // above the blue isolated vertex there may be one blue subcurve, but
      // the subcruve above it must be red in this case. This is why it is
      // sufficient to go at most two steps up.
      // There is nothing above the vertex - use the current red top face.
      Subcurve* sc_above = *iter;
      if (sc_above == NULL) {
        red_face = m_overlay_helper.red_top_face();
      }
      else {
        if (sc_above->color() != Traits_2::BLUE) {
          red_face = sc_above->red_halfedge_handle()->face();
        }
        else {
          sc_above = sc_above->subcurve_above();
          red_face = (sc_above != NULL) ?
            sc_above->red_halfedge_handle()->face() :
            m_overlay_helper.red_top_face();
        }
      }
    }

    // Use the overlay traits to update the newly created isolated vertex.
    m_overlay_traits->create_vertex(red_face, *blue_vertex, new_v);
  }
  else {
    // We have an isolated red vertex inside a blue face.
    Face_handle_blue blue_face;

    // Obtain the blue face containing the isolated vertex.
    if (iter == this->status_line_end()) {
      // There is nothing above the vertex - use the current blue top face.
      blue_face = m_overlay_helper.blue_top_face();
    }
    else {
      // Go up the status line until we find a blue halfedge.
      // Note that the subcurve above is always of an opposite color. However,
      // above the red isolated vertex there may be one red subcurve, but
      // the subcruve above it must be blue in this case. This is why it is
      // sufficient to go at most two steps up.
      // If we do not find a blue halfedge, we use the current red top face.
      Subcurve* sc_above = *iter;
      if (sc_above == NULL) {
        blue_face = m_overlay_helper.blue_top_face();
      }
      else {
        if (sc_above->color() != Traits_2::RED) {
          blue_face = sc_above->blue_halfedge_handle()->face();
        }
        else {
          sc_above = sc_above->subcurve_above();
          blue_face = (sc_above != NULL) ?
            sc_above->blue_halfedge_handle()->face() :
            m_overlay_helper.blue_top_face();
        }
      }
    }

    // Use the overlay traits to update the newly created isolated vertex.
    m_overlay_traits->create_vertex(*red_vertex, blue_face, new_v);
  }

  return new_v;
}

//-----------------------------------------------------------------------------
// Map a newly created halfedge in the result arrangement to its originator
// red and blue halfedges.
//
template <typename OvlHlpr, typename OvlTr>
void Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
_map_halfedge_and_twin(Halfedge_handle he,
                       Halfedge_handle_red red_he,
                       Halfedge_handle_blue blue_he)
{
  if (he->direction() == ARR_LEFT_TO_RIGHT) he = he->twin();

  // Obtain the twin red and blue halfedges (if they are valid). Note that
  // the original halfedges are always directed from right to left.
  Halfedge_handle_red     red_he_twin;
  Halfedge_handle_blue    blue_he_twin;

  if (red_he != Halfedge_handle_red()) red_he_twin = red_he->twin();
  if (blue_he != Halfedge_handle_blue()) blue_he_twin = blue_he->twin();

  // Create the halfedge-information pairs and store them in the map.
  Halfedge_info he_info = std::make_pair(red_he, blue_he);
  Halfedge_info he_twin_info = std::make_pair(red_he_twin, blue_he_twin);

  m_halfedges_map[he] = he_info;
  m_halfedges_map[he->twin()] = he_twin_info;
}

//-----------------------------------------------------------------------------
// Update the boundary vertices map.
//
template <typename OvlHlpr, typename OvlTr>
void Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
_map_boundary_vertices(Event* event, Vertex_handle v, boost::mpl::bool_<true>)
{
  // Update the red and blue object if the last event on sc is on the boundary.
  if ((event->parameter_space_in_x() != ARR_INTERIOR) ||
      (event->parameter_space_in_y() != ARR_INTERIOR))
  {
    const Point_2& pt = event->point();
    typename Vertex_map::iterator it = m_vertices_map.find(v);
    if (it == m_vertices_map.end()) return;
    Handle_info& info = it->second;
    const Cell_handle_blue* blue_handle_p = pt.blue_cell_handle();
    if (blue_handle_p) info.second = *blue_handle_p;
    const Cell_handle_red* red_handle_p = pt.red_cell_handle();
    if (red_handle_p) info.first = *red_handle_p;

    if (!boost::get<Face_handle_red>(&(info.first)) &&
        !boost::get<Face_handle_blue>(&(info.second)))
    {
      // If both, the red and blue, variants do not represent face handles,
      // they must represt either vertex or edge handles. In this case it is
      // safe to apply the call to the overlay traits and erase the record,
      // so that it is not applied again for this vertex.
      const Cell_handle_red& red_handle = info.first;
      const Cell_handle_blue& blue_handle = info.second;
      Vertex_handle v = (*it).first;
      Create_vertex_visitor visitor(m_overlay_traits, v);
      boost::apply_visitor(visitor, red_handle, blue_handle);
      m_vertices_map.erase(it);
    }
  }
}

//-----------------------------------------------------------------------------
// Update the boundary vertices map.
//
template <typename OvlHlpr, typename OvlTr>
void Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
_map_boundary_vertices(Event* /* event */, Vertex_handle /* v */, boost::mpl::bool_<false>)
{}

/* Notify the overlay traits about a newly created vertex.
 * If the new vertex is on the boundary and the red object (or blue object) is
 * empty, the red object (or blue object) might have not been updated yet. In
 * this case, we postpone the notification for all (contracted and identified)
 * boundary side vertices to the end of the sweep.
 */
template <typename OvlHlpr, typename OvlTr>
void Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
_create_vertex(Event* event,
               Vertex_handle new_v,
               Subcurve* sc,
               boost::mpl::bool_<true>)
{
  const Point_2& pt = event->point();
  const Cell_handle_red* red_handle = pt.red_cell_handle();
  const Cell_handle_blue* blue_handle = pt.blue_cell_handle();

  // If the vertex is on the boundary, postpone the notification, but
  // update the red and objects in case they are empty.
  if ((event->parameter_space_in_x() != ARR_INTERIOR) ||
      (event->parameter_space_in_y() != ARR_INTERIOR))
  {
    if (!red_handle) {
      CGAL_assertion(blue_handle != NULL);
      // Obtain the red face by looking for a subcurve above.
      const Subcurve* sc_above = sc->subcurve_above();
      Face_handle_red red_f = (sc_above != NULL) ?
        sc_above->red_halfedge_handle()->face() : sc->red_top_face();
      Handle_info info = std::make_pair(Cell_handle_red(red_f), *blue_handle);
      m_vertices_map[new_v] = info;
      return;
    }
    if (!blue_handle) {
      CGAL_assertion(red_handle != NULL);
      // Obtain the blue face by looking for a subcurve above.
      const Subcurve* sc_above = sc->subcurve_above();
      Face_handle_blue blue_f = (sc_above != NULL) ?
        sc_above->blue_halfedge_handle()->face() : sc->blue_top_face();
      Handle_info info = std::make_pair(*red_handle, Cell_handle_blue(blue_f));
      m_vertices_map[new_v] = info;
      return;
    }

    Handle_info info = std::make_pair(*red_handle, *blue_handle);
    m_vertices_map[new_v] = info;
    return;
  }

  _create_vertex(event, new_v, sc, boost::mpl::bool_<false>());
}

/* Notify the overlay traits about a newly created vertex. */
template <typename OvlHlpr, typename OvlTr>
void Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::
_create_vertex(Event* event,
               Vertex_handle new_v,
               Subcurve* sc,
               boost::mpl::bool_<false>)
{
  const Point_2& pt = event->point();
  const Cell_handle_red* red_handle = pt.red_cell_handle();
  const Cell_handle_blue* blue_handle = pt.blue_cell_handle();
  if (!red_handle) {
    // A blue vertex is located inside a red face. Obtain the red face
    // by looking for a subcurve above.
    const Subcurve* sc_above = sc->subcurve_above();
    Face_handle_red red_f = (sc_above != NULL) ?
      sc_above->red_halfedge_handle()->face() : sc->red_top_face();

    CGAL_assertion(blue_handle != NULL);
    const Vertex_handle_blue& blue_v =
      boost::get<Vertex_handle_blue>(*blue_handle);
    m_overlay_traits->create_vertex(red_f, blue_v, new_v);
    return;
  }

  if (!blue_handle) {
    // A red vertex is located inside a blue face. Obtain the blue face
    // by looking for a subcurve above.
    const Subcurve* sc_above = sc->subcurve_above();
    Face_handle_blue blue_f = (sc_above != NULL) ?
      sc_above->blue_halfedge_handle()->face() : sc->blue_top_face();

    CGAL_assertion(red_handle != NULL);
    const Vertex_handle_red& red_v =
      boost::get<Vertex_handle_red>(*red_handle);
    m_overlay_traits->create_vertex(red_v, blue_f, new_v);
    return;
  }

  Create_vertex_visitor visitor(m_overlay_traits, new_v);
  boost::apply_visitor(visitor, *red_handle, *blue_handle);
}

//-----------------------------------------------------------------------------
// Update a newly created result edge using the overlay traits.
//
template <typename OvlHlpr, typename OvlTr>
void
Arr_overlay_sl_visitor<OvlHlpr, OvlTr>::_create_edge(Subcurve* sc,
                                                     Halfedge_handle new_he)
{
  // Note that the "red" and "blue" halfedges are always directed from right
  // to left, so we make sure the overlaid halfedge is also directed from
  // right to left.
  if (new_he->direction() != ARR_RIGHT_TO_LEFT) new_he = new_he->twin();

  // Examine the various cases for the creation of a new edge.
  if (sc->color() == Traits_2::RB_OVERLAP) {
    // The new edge represents an overlap between a red halfedge and a blue
    // halfedge.
    Halfedge_handle_red   red_he = sc->red_halfedge_handle();
    Halfedge_handle_blue  blue_he = sc->blue_halfedge_handle();

    m_overlay_traits->create_edge(red_he, blue_he, new_he);
  }
  else if (sc->color() == Traits_2::RED) {
    // We have a red edge on a blue face.
    Halfedge_handle_red   red_he = sc->red_halfedge_handle();
    Subcurve* sc_above = sc->subcurve_above();
    Face_handle_blue blue_f = (sc_above != NULL) ?
      sc_above->blue_halfedge_handle()->face() : sc->blue_top_face();
    m_overlay_traits->create_edge(red_he, blue_f, new_he);
  }
  else {
    CGAL_assertion(sc->color() == Traits_2::BLUE);

    // We have a blue edge on a red face.
    Halfedge_handle_blue  blue_he = sc->blue_halfedge_handle();
    Subcurve* sc_above = sc->subcurve_above();
    Face_handle_red red_f = (sc_above != NULL) ?
      sc_above->red_halfedge_handle()->face() : sc->red_top_face();
    m_overlay_traits->create_edge(red_f, blue_he, new_he);
  }
}

} //namespace CGAL

#endif
