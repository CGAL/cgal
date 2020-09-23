// Copyright (c) 2006,2007,2008,2009,2010,2011,2012,2013 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_ARR_CONSTRUCTION_SS_VISITOR_H
#define CGAL_ARR_CONSTRUCTION_SS_VISITOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#ifndef CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
#define CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE 0
#endif

/*! \file
 *
 * Definition of the Arr_construction_ss_visitor class-template.
 */

#include <vector>

#include <CGAL/Arr_accessor.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Surface_sweep_2/Default_visitor_base.h>
#include <CGAL/Default.h>

namespace CGAL {

namespace Ss2 = Surface_sweep_2;

/*! \struct Integer_hash_function
 * An auxiliary hash functor for integers.
 */
struct Integer_hash_function {
  typedef std::size_t result_type;

  std::size_t operator()(unsigned int i) const { return i; }
};

/*! \class Arr_construction_ss_visitor
 * A sweep-line visitor for constructing an arrangement embedded on a surface.
 */
template <typename Helper_, typename Visitor_ = Default>
class Arr_construction_ss_visitor :
  public Ss2::Default_visitor_base<typename Helper_::Geometry_traits_2,
                                   typename Helper_::Event,
                                   typename Helper_::Subcurve,
                                   typename Helper_::Allocator,
                                   typename Default::Get<
                                     Visitor_,
                                     Arr_construction_ss_visitor<
                                       Helper_, Visitor_> >::type>
{
public:
  typedef Helper_                                       Helper;

  typedef typename Helper::Geometry_traits_2            Geometry_traits_2;
  typedef typename Helper::Event                        Event;
  typedef typename Helper::Subcurve                     Subcurve;
  typedef typename Helper::Allocator                    Allocator;

private:
  typedef Geometry_traits_2                             Gt2;
  typedef Arr_construction_ss_visitor<Helper, Visitor_> Self;
  typedef typename Default::Get<Visitor_, Self>::type   Visitor;
  typedef Ss2::Default_visitor_base<Gt2, Event, Subcurve, Allocator, Visitor>
                                                        Base;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

protected:
  typedef typename Helper::Arrangement_2                Arrangement_2;
  typedef typename Arrangement_2::Topology_traits       Topology_traits;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;

  typedef typename Base::Event_subcurve_iterator        Event_subcurve_iterator;
  typedef typename Base::Event_subcurve_reverse_iterator
    Event_subcurve_reverse_iterator;
  typedef typename Subcurve::Status_line_iterator       Status_line_iterator;

  typedef typename Helper::Indices_list                 Indices_list;
  typedef typename Helper::Halfedge_indices_map         Halfedge_indices_map;
  typedef Unique_hash_map<unsigned int, Vertex_handle, Integer_hash_function>
                                                        Iso_vertices_map;

protected:
  Helper m_helper;                          // The helper class.

  Arrangement_2* m_arr;                     // The arrangement we construct.
  Topology_traits* m_top_traits;            // The topology-traits class.
  Arr_accessor<Arrangement_2> m_arr_access; // An arrangement accessor.

  unsigned int m_sc_counter;                // Counter for subcurves that may
                                            // represent a hole (the upper
                                            // subcurves that emarge from event
                                            // points with only right curves).

  std::vector<Halfedge_handle> m_sc_he_table; // A table that maps a subcurve
                                            // index to its halfedge handle,
                                            // directed from right to left.

  Iso_vertices_map m_iso_verts_map;         // Maps an index to the isolated
                                            // vertex.

  Halfedge_indices_map m_he_indices_table;  // Maps each halfdge to the
                                            // indices of subcurves that
                                            // lies below it.

  const Vertex_handle m_invalid_vertex;     // An invalid vertex handle.

public:
  /*! Constructor. */
  Arr_construction_ss_visitor(Arrangement_2* arr) :
    m_helper(arr),
    m_arr(arr),
    m_top_traits(arr->topology_traits()),
    m_arr_access(*arr),
    m_sc_counter(0),
    m_sc_he_table(1),
    m_invalid_vertex()
  { m_helper.set_halfedge_indices_map(m_he_indices_table); }

  /*! Destructor. */
  virtual ~Arr_construction_ss_visitor() {}

  /// \name Sweep-line notifications.
  //@{

  /* A notification issued before the sweep process starts. */
  inline void before_sweep();

  /*!
   * A notification invoked before the sweep-line starts handling the given
   * event.
   */
  inline void before_handle_event(Event* event);

  /*!
   * A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  bool after_handle_event(Event* event, Status_line_iterator iter, bool flag);

  /*! A notification invoked when a new subcurve is created. */
  void add_subcurve(const X_monotone_curve_2& cv, Subcurve* sc);
  //@}

  /// \name Insertion functions.
  //@{

  /*!
   * Insert the given subcurve in the interior of a face.
   * \param cv The geometric subcurve.
   * \param sc The sweep-line subcurve information.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle
  insert_in_face_interior(const X_monotone_curve_2& cv, Subcurve* sc);

  /*!
   * Insert the given subcurve given its left end-vertex.
   * \param cv The geometric entity.
   * \param prev The predecessor halfedge around the left vertex.
   * \param sc The sweep-line subcurve information.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle
  insert_from_left_vertex(const X_monotone_curve_2& cv, Halfedge_handle he,
                          Subcurve* sc);

  /*!
   * Insert the given subcurve given its right end-vertex.
   * \param cv The geometric entity.
   * \param prev The predecessor halfedge around the right vertex.
   * \param sc The sweep-line subcurve information.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle
  insert_from_right_vertex(const X_monotone_curve_2& cv,
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

  /*!
   * Relocate holes and isolated vertices inside a newly created face f2,
   * that was split from f1 after the insertion of a new edge.
   * \param he The halfedge that caused the face split. Its incident face is
   *           the new face f2, and the incident face of its twin is f1.
   */
  void relocate_in_new_face(Halfedge_handle he);
  //@}

  /*! Get the last event associated with the given subcurve. */
  Event* last_event_on_subcurve(Subcurve* sc) { return sc->last_event(); }

private:
  /// \name Auxiliary functions.
  //@{

  /*!
   * Cast a Traits::Point_2 object into an Arrangement_2::Point_2 object.
   * These two types may not be the same when the addition visitor inherits
   * from this base class.
   */
  inline const typename Arrangement_2::Point_2& _point(const Point_2& p) const
  { return (static_cast<const typename Arrangement_2::Point_2&>(p)); }

  /*!
   * Cast a Traits::X_monotone_curve_2 object into an
   * Arrangement_2::X_monotone_curve_2 object.
   * These two types may not be the same when the addition visitor inherits
   * from this base class.
   */
  inline const typename Arrangement_2::X_monotone_curve_2&
  _curve(const X_monotone_curve_2& cv) const
  {
    return (static_cast<const typename Arrangement_2::X_monotone_curve_2&>(cv));
  }

  /*! Map the given subcurve index to the given halfedge handle. */
  void _map_new_halfedge(unsigned int i, Halfedge_handle he);
  //@}
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
// Notifies the helper that the sweep process now starts.
template <typename Hlpr, typename Vis>
void Arr_construction_ss_visitor<Hlpr, Vis>::before_sweep()
{ m_helper.before_sweep(); }

//-----------------------------------------------------------------------------
// A notification invoked before the sweep-line starts handling the given
// event.
//
template <class Hlpr, typename Vis>
void Arr_construction_ss_visitor<Hlpr, Vis>::before_handle_event(Event* event)
{
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << "CGAL_CSLV before_handle_event" << std::endl;
#endif
  // We just have to notify the helper class on the event.
  m_helper.before_handle_event(event);
}

//-----------------------------------------------------------------------------
// A notification invoked after the sweep-line finishes handling the given
// event.
//
template <typename Hlpr, typename Vis>
bool Arr_construction_ss_visitor<Hlpr, Vis>::
after_handle_event(Event* event, Status_line_iterator iter, bool /* flag */)
{
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << std::endl << "CGAL_CSLV after_handle_event: " << event
            << std::endl;
#endif

  // Check if the event represents an isolated vertex.
  if (!event->has_left_curves() && !event->has_right_curves()) {
    // There are no incident subcurves, so this event is an isolated vertex.
    // We map the current index to this vertex, and add this index to the
    // indices list of the curve the vertex "sees" from below.
    Vertex_handle v = insert_isolated_vertex(event->point(), iter);

    ++m_sc_counter;
    m_iso_verts_map[m_sc_counter] = v;
    _map_new_halfedge(m_sc_counter, Halfedge_handle());

    if (iter != this->status_line_end()) {
      // The isolated vertex "sees" the subcurve of the given position from
      // below.
      Subcurve* sc_above = *iter;
      sc_above->add_halfedge_index(m_sc_counter);
    }
    else {
      // The vertex is not located below any valid curve, so we use the helper
      // class to mark that this index should belong to the current top face.
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
      std::cout << "CGAL_CSLV adding a " << m_sc_counter << std::endl;
#endif
      m_helper.add_subcurve_in_top_face(m_sc_counter);
    }

    // The event can now be deallocated.
    return true;
  }

  // TODO EBEB 2012-10-16 compile only when non-oblivious
  if (event->parameter_space_in_x() == CGAL::ARR_LEFT_BOUNDARY) {
    if (!this->is_status_line_empty()) {
      Status_line_iterator prev = iter;
      for (size_t i = 0; i < event->number_of_right_curves(); ++i) --prev;
      // move items from top face to last inserted curve
      Indices_list& list_ref = (*prev)->halfedge_indices_list();
      list_ref.clear();
      list_ref.splice(list_ref.end(), m_helper.halfedge_indices_list());

#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
      typename Indices_list::const_iterator it;
      for (it = list_ref.begin(); it != list_ref.end(); ++it)
        std::cout << "moved " << *it << " from top to below" << std::endl;
#endif
    }
  }

  // Check if the event has only incident subcurves from its right.
  if (!event->has_left_curves()) {
    CGAL_assertion(event->has_right_curves());

    // In case of a finite event that has no incident left curves, it is
    // associated with a point that may be the leftmost one in a hole.
    // We give index to the topmost subcurve from the right, and add this
    // vertex indices list of the curve the event "sees" from below.
    ++m_sc_counter;
    (*(event->right_curves_rbegin()))->set_index(m_sc_counter);
    if (iter != this->status_line_end()) {
      // The vertex "sees" the subcurve of the given position from below.
      Subcurve* sc_above = *iter;
      sc_above->add_halfedge_index(m_sc_counter);
    }
    else {
      // The vertex is not located below any valid curve, so we use the helper
      // class to mark that this index should belong to the current top face.
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
      std::cout << "CGAL_CSLV adding b " << m_sc_counter << std::endl;
#endif
      m_helper.add_subcurve_in_top_face(m_sc_counter);
    }
  }

  // Set the last event of all left subcurve (thus, this event corresponds
  // to their right endpoint).
  Event_subcurve_iterator  left_it;
  for (left_it = event->left_curves_begin();
       left_it != event->left_curves_end(); ++left_it)
    (*left_it)->set_last_event(event);

  // In case there are no right subcurves, the event can be deallocated.
  if (event->number_of_right_curves() == 0) {
    // Inform the helper class that the event will soon be deallocated.
    m_helper.before_deallocate_event(event);
    return true;
  }

  // Mark that all right subcurves incident to the current event are not
  // in the arrangement yet.
  event->init_subcurve_in_arrangement_flags(event->number_of_right_curves());

  // Set the last event of all right subcurve (thus, this event corresponds
  // to their left endpoint).
  Event_subcurve_iterator  right_it;
  for (right_it = event->right_curves_begin();
       right_it != event->right_curves_end(); ++right_it)
    (*right_it)->set_last_event(event);

  // Mark that the event cannot be deallocated just yet.
  return false;
}

//-----------------------------------------------------------------------------
// A notification invoked when a new subcurve is created.
//
template <typename Hlpr, typename Vis>
void Arr_construction_ss_visitor<Hlpr, Vis>::
add_subcurve(const X_monotone_curve_2& cv, Subcurve* sc)
{
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << std::endl << "CGAL_CSLV add_subcurve: " << cv << std::endl;
#endif

  // Obtain all information to perform the insertion of the subcurve into
  // the arrangement.
  Event* last_event = last_event_on_subcurve(sc);
  Halfedge_handle res;
  Halfedge_handle he_right = this->current_event()->halfedge_handle();
  Halfedge_handle he_left = last_event->halfedge_handle();
  const int jump = last_event->compute_halfedge_jump_count(sc);

  const Halfedge_handle invalid_he;

#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  if (last_event->is_closed())
    std::cout << "CGAL_CSLG lastevent: " << last_event->point() << std::endl;
  if (he_left != invalid_he) {
    std::cout << "he_left    : " << &(*he_left) << std::endl;
    if (!he_left->is_fictitious())
      std::cout << "he_leftcv  : " << he_left->curve() << std::endl;
    else std::cout << "he_left    : fictitious" << std::endl;
    std::cout << "he_leftdir : " << he_left->direction() << std::endl;
    std::cout << "he_leftfac : " << &(*he_left->face()) << std::endl;
  }
  else std::cout << "he_left    : invalid" << std::endl;
  if (he_right != invalid_he) {
    std::cout << "he_right   : " << &(*he_right) << std::endl;
    if (!he_right->is_fictitious())
      std::cout << "he_rightcv : " << he_right->curve() << std::endl;
    else std::cout << "he_right   : fictitious" << std::endl;
    std::cout << "he_rightdir: " << he_right->direction() << std::endl;
    std::cout << "he_rightfac: " << &(*he_right->face()) << std::endl;
  } else std::cout << "he_right   : invalid" << std::endl;
#endif

  // Check whether the previous event on the curve is not in the arrangement
  // yet.
  if (he_left == invalid_he) {
    Vertex_handle v_left = last_event->vertex_handle();
    // Check whether the vertex to be associated with the left end of
    // the curve has already been created.
    if ((v_left != m_invalid_vertex) && (v_left->degree() > 0)) {
      // The left vertex v is a boundary vertex which already has some
      // incident halfedges. We look for the predecessor halfedge and
      // insert the curve between two existing vertices.
      Arr_parameter_space bx = last_event->parameter_space_in_x();
      Arr_parameter_space by = last_event->parameter_space_in_y();
      CGAL_assertion((bx != ARR_INTERIOR) || (by != ARR_INTERIOR));
      he_left = Halfedge_handle
        (m_top_traits->locate_around_boundary_vertex(&(*v_left), _curve(cv),
                                                     ARR_MIN_END, bx, by));
    }
  }
  else {
    // The previous event on the curve is already in the arrangement.
    // Therefore, we use it to insert the subcurve.
    // First, we skip some halfedges around the left vertex to get the true
    // predecessor halfedge for the insertion.
    for (int i = 0; i < jump; i++) he_left = (he_left->next())->twin();

#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
    if (jump != 0) {
      std::cout << "CGAL_CSLV JUMP: " << jump << std::endl;
      if (!he_left->is_fictitious())
        std::cout << "he_leftcv  : " << he_left->curve() << std::endl;
      else std::cout << "he_left    : fictitious" << std::endl;
      std::cout << "he_leftdir : " << he_left->direction() << std::endl;
      std::cout << "he_leftfac : " << &(*he_left->face()) << std::endl;
    }
#endif
  }

  // Check whether the current event is already in the arrangement
  if (he_right == invalid_he) {
    // We do not have handles for any of the curve end, so we insert it in
    // the interior of a face.
    Event* curr_event = this->current_event();
    Vertex_handle v_right = curr_event->vertex_handle();
    if ((v_right != m_invalid_vertex) && (v_right->degree() > 0)) {
      // The left vertex v is a boundary vertex which already has some
      // incident halfedges. We look for the predecessor halfedge and
      // insert the curve from this right vertex.
      Arr_parameter_space bx = curr_event->parameter_space_in_x();
      Arr_parameter_space by = curr_event->parameter_space_in_y();
      CGAL_assertion((bx != ARR_INTERIOR) || (by != ARR_INTERIOR));
      he_right = Halfedge_handle
        (m_top_traits->locate_around_boundary_vertex(&(*v_right), _curve(cv),
                                                     ARR_MAX_END, bx, by));
    }
  }

  // Check whether the verteices to be associated with the left end and
  // the right end of the curve have already been created.
  bool dummy;
  res = (he_left != invalid_he) ?
    ((he_right != invalid_he) ?
     this->insert_at_vertices(cv, he_right, he_left, sc, dummy) :
     this->insert_from_left_vertex(cv, he_left, sc)) :
    ((he_right != invalid_he) ?
     this->insert_from_right_vertex(cv, he_right, sc) :
     this->insert_in_face_interior(cv, sc));

#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << "CGAL_CSLV res: " << &(*res) << " with face " << &(*res->face())
            << " is " << res->direction() << std::endl;
  std::cout << "CGAL_CSLV twi: " << &(*res->twin()) << " with face "
            << &(*res->twin()->face()) << " is " << res->twin()->direction()
            << std::endl;
#endif

  // Make sure that res is a halfedge that is always directed from left to
  // right (thus its twin is directed from right to left).
  if (res->direction() != ARR_LEFT_TO_RIGHT) {
    res = res->twin();
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
    std::cout << "CGAL_CSLV twined!" << std::endl;
#endif
  }

  // Update the last event with the inserted halfegde (if necessary)
  // and check if we have to update the auxiliary information on the location
  // of holes.
  if ((last_event->number_of_left_curves() == 0) &&
      last_event->is_curve_largest((Subcurve*)sc))
  {
    if (last_event->vertex_handle() == m_invalid_vertex)
      last_event->set_halfedge_handle(res->twin());

    // If sc has valid index, insert its index to m_sc_he_table.
    if (sc->has_valid_index()) {
      CGAL_assertion(res->twin()->direction() == ARR_RIGHT_TO_LEFT);
      _map_new_halfedge(sc->index(), res->twin());
    }
  }

  // Update the halfedge handle associated with the current event.
  if (this->current_event()->vertex_handle() == m_invalid_vertex)
      this->current_event()->set_halfedge_handle(res);

  // In case the event has no more right subcurves associated with it, we can
  // deallocate it. Note that we inform the helper class before deallocating
  // the event.
  if (((Event*) sc->right_event())==this->current_event() &&
      last_event->dec_right_curves_counter() == 0)
  {
    m_helper.before_deallocate_event(last_event);
    this->deallocate_event(last_event);
  }

  // Clear the list of indices of the subcurve.
  sc->clear_halfedge_indices();
}

//-----------------------------------------------------------------------------
// Insert the given subcurve in the interior of an arrangement face.
//
template <typename Hlpr, typename Vis>
typename Arr_construction_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_construction_ss_visitor<Hlpr, Vis>::
insert_in_face_interior(const X_monotone_curve_2& cv, Subcurve* sc)
{
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << "CGAL_CSLV insert_in_face_interior\ncurve: " << cv << std::endl;
#endif

  Event* last_event = last_event_on_subcurve(sc);
  Vertex_handle  v1 = last_event->vertex_handle();
  CGAL_assertion((v1 == m_invalid_vertex) || (v1->degree() == 0));

  if (v1 == m_invalid_vertex)
    // Create the vertex to be associated with the left end of the curve.
    v1 = m_arr_access.create_vertex(_point(last_event->point()));

  Event* curr_event = this->current_event();
  Vertex_handle v2 = curr_event->vertex_handle();
  CGAL_assertion((v2 == m_invalid_vertex) || (v2->degree() == 0));

  if (v2 == m_invalid_vertex)
    // Create the vertex to be associated with the right end of the curve.
    v2 = m_arr_access.create_vertex(_point(curr_event->point()));

  // Perform the insertion between the two (currently isolated) vertices in
  // the interior of the current top face, as given by the helper class.
  Halfedge_handle res =
    m_arr_access.insert_in_face_interior_ex(m_helper.top_face(), _curve(cv),
                                            ARR_LEFT_TO_RIGHT, v1, v2);

  // Map the new halfedge to the indices list of all subcurves that lie
  // below it.
  if (sc->has_halfedge_indices()) {
    CGAL_assertion(res->twin()->direction() == ARR_RIGHT_TO_LEFT);
    Indices_list& list_ref = m_he_indices_table[res->twin()];
    list_ref.clear();
    list_ref.splice(list_ref.end(), sc->halfedge_indices_list());
  }

  // Notify the helper on the creation of the new halfedge.
  m_helper.add_subcurve(res, sc);

  return res;
}

//-----------------------------------------------------------------------------
// Insert the given subcurve using its two end-vertices.
//
template <typename Hlpr, typename Vis>
typename Arr_construction_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_construction_ss_visitor<Hlpr, Vis>::
insert_at_vertices(const X_monotone_curve_2& cv,
                   Halfedge_handle prev1,
                   Halfedge_handle prev2,
                   Subcurve* sc,
                   bool& new_face_created)
{
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << "CGAL_CSLV insert_at_vertices:\ncurve:" << cv << std::endl;
  if (!prev1->is_fictitious())
    std::cout << "prev1cv   : " << prev1->curve() << std::endl;
  else std::cout << "prev1     : fictitious" << std::endl;
  std::cout << "prev1dir  : " << prev1->direction() << std::endl;
  std::cout << "prev1fac  : " << &(*prev1->face()) << std::endl;
  if (!prev2->is_fictitious())
    std::cout << "prev2cv   : " << prev2->curve() << std::endl;
  else std::cout << "prev2     : fictitious" << std::endl;
  std::cout << "prev2dir  : " << prev2->direction() << std::endl;
  std::cout << "prev2fac  : " << &(*prev2->face()) << std::endl;
#endif

  // Use the helper class to determine whether the order of predecessor
  // halfedges should be swaped, to that the edge directed from prev1->target()
  // to prev2->target() is incident to the new face (in case a new face is
  // created).
  Halfedge_handle res;
  const bool swap_preds = m_helper.swap_predecessors(this->current_event());

  // Comment: In some topologies swap_preds is always false,
  //          thus we use 'false' to disallow swapping
  //          such that the compiler can optimize awy the
  //          computation of signs and local/global minima
  //          that now takes place inside _insert_at_vertices
  // TODO EBEB 2012-08-06 check whether signs are not needed
  // it seems that swap_pred is either
  //   false
  // if oblivious or
  //   event->parameter_space_in_x() == CGAL::ARR_INTERIOR &&
  //   event->parameter_space_in_y() == CGAL::ARR_TOP_BOUNDARY
  // if not oblivious! But I have the feeling that signs are needed!

  bool check_swapped_predecessors = true;
#if 1
  res = (swap_preds) ?
    // if swapping prev2 will become outer of new split face (it it exists)
    m_arr_access.insert_at_vertices_ex(prev2, _curve(cv), ARR_LEFT_TO_RIGHT,
                                       prev1->next(), new_face_created,
                                       check_swapped_predecessors, false) :
    // usually prev1 is outer of new split face (it it exists)
    // order is determined by top-traits helper!
    // "false" disallows swapping of prev1/preve2! ...
    m_arr_access.insert_at_vertices_ex(prev1, _curve(cv), ARR_RIGHT_TO_LEFT,
                                       prev2->next(), new_face_created,
                                       check_swapped_predecessors, false);

  // ... thus the value should now have changed
  CGAL_assertion(!check_swapped_predecessors);
#else
  res = m_arr_access.insert_at_vertices_ex(prev1, _curve(cv), ARR_RIGHT_TO_LEFT,
                                           prev2->next(), new_face_created,
                                           check_swapped_predecessors);
  if (check_swapped_predecessors) res = res->twin();
#endif

  // Map the new halfedge to the indices list of all subcurves that lie
  // below it.
  if (sc->has_halfedge_indices()) {
    Halfedge_handle he = res;
    if (swap_preds) he = he->twin();
    CGAL_assertion(he->direction() == ARR_RIGHT_TO_LEFT);
    Indices_list& list_ref = m_he_indices_table[he];
    list_ref.clear();
    list_ref.splice(list_ref.end(), sc->halfedge_indices_list());
  }

  // Notify the helper on the creation of the new halfedge.
  // Note that we do this before trying to relocate holes in the new
  // face (if one is created).
  m_helper.add_subcurve(res, sc);

  if (new_face_created) {
    // TODO EBEB 2012-10-17 needs a tests for outer-outer insertion
    // In case a new face has been created (pointed by the new halfedge
    // we obtained), we have to examine the holes and isolated vertices
    // in the existing face (pointed by the twin halfedge) and relocate
    // the relevant features in the new face.
    CGAL_assertion(res->face() != res->twin()->face());
    this->relocate_in_new_face(res);
  }

  return res;
}

//-----------------------------------------------------------------------------
// Insert the given subcurve from a vertex that corresponds to its right end.
//
template <typename Hlpr, typename Vis>
typename Arr_construction_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_construction_ss_visitor<Hlpr, Vis>::
insert_from_right_vertex(const X_monotone_curve_2& cv,
                         Halfedge_handle prev,
                         Subcurve* sc)
{
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << "CGAL_CSLV insert_from_right_vertex:\ncurve:" << cv << std::endl;
  if (!prev->is_fictitious()) {
    std::cout << "prevcv    : " << prev->curve() << std::endl;
  } else {
    std::cout << "prev      : fictitious" << std::endl;
  }
  std::cout << "prevdir   : " << prev->direction() << std::endl;
  std::cout << "prevfac   : " << &(*prev->face()) << std::endl;
#endif

  Event* last_event = last_event_on_subcurve(sc);
  Vertex_handle v = last_event->vertex_handle();
  CGAL_assertion((v == m_invalid_vertex) || (v->degree() == 0));

  // Create the vertex to be associated with the left end of the curve.
  if (v == m_invalid_vertex)
    v = m_arr_access.create_vertex(_point(last_event->point()));

  // Insert the curve given its left vertex and the predecessor around the
  // right vertex.
  Halfedge_handle res =
    m_arr_access.insert_from_vertex_ex(prev, _curve(cv), ARR_RIGHT_TO_LEFT, v);

  // Map the new halfedge to the indices list of all subcurves that lie
  // below it.
  if (sc->has_halfedge_indices()) {
    CGAL_assertion(res->direction() == ARR_RIGHT_TO_LEFT);
    Indices_list& list_ref = m_he_indices_table[res];
    list_ref.clear();
    list_ref.splice(list_ref.end(), sc->halfedge_indices_list());
  }

  // Notify the helper on the creation of the new halfedge.
  m_helper.add_subcurve(res, sc);
  return res;
}

//-----------------------------------------------------------------------------
// Insert the given subcurve from a vertex that corresponds to its left end.
//
template <typename Hlpr, typename Vis>
typename Arr_construction_ss_visitor<Hlpr, Vis>::Halfedge_handle
Arr_construction_ss_visitor<Hlpr, Vis>::
insert_from_left_vertex(const X_monotone_curve_2& cv,
                        Halfedge_handle prev,
                        Subcurve* sc)
{
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << "CGAL_CSLV insert_from_left_vertex:\ncurve:" << cv << std::endl;
  if (!prev->is_fictitious())
    std::cout << "prevcv    : " << prev->curve() << std::endl;
  else std::cout << "prev      : fictitious" << std::endl;
  std::cout << "prevdir   : " << prev->direction() << std::endl;
  std::cout << "prevfac   : " << &(*prev->face()) << std::endl;
#endif

  Event* curr_event = this->current_event();
  Vertex_handle v = curr_event->vertex_handle();
  CGAL_assertion((v == m_invalid_vertex) || (v->degree() == 0));

  // Create the vertex to be associated with the right end of the curve.
  if (v == m_invalid_vertex)
    v = m_arr_access.create_vertex(_point(curr_event->point()));

  // Insert the curve given its right vertex and the predecessor around the
  // left vertex.
  Halfedge_handle  res =
    m_arr_access.insert_from_vertex_ex(prev, _curve(cv), ARR_LEFT_TO_RIGHT, v);

  // Map the new halfedge to the indices list of all subcurves that lie
  // below it.
  if (sc->has_halfedge_indices()) {
    CGAL_assertion(res->twin()->direction() == ARR_RIGHT_TO_LEFT);
    Indices_list& list_ref = m_he_indices_table[res->twin()];
    list_ref.clear();
    list_ref.splice(list_ref.end(), sc->halfedge_indices_list());
  }

  // Notify the helper on the creation of the new halfedge.
  m_helper.add_subcurve(res, sc);
  return res;
}

//-----------------------------------------------------------------------------
// Insert an isolated vertex into the arrangement.
//
template <typename Hlpr, typename Vis>
typename Arr_construction_ss_visitor<Hlpr, Vis>::Vertex_handle
Arr_construction_ss_visitor<Hlpr, Vis>::
insert_isolated_vertex(const Point_2& pt, Status_line_iterator /* iter */)
{
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << "CGAL_CSLV insert_isolated_vertex:\npoint:" << pt << std::endl;
#endif

  // Insert the isolated vertex in the interior of the current top face, as
  // given by the helper class.
  return m_arr->insert_in_face_interior(_point(pt), m_helper.top_face());
}

//-----------------------------------------------------------------------------
// Reloacte holes and isolated vertices inside a newly created face.
//
template <typename Hlpr, typename Vis>
void Arr_construction_ss_visitor<Hlpr, Vis>::
relocate_in_new_face(Halfedge_handle he)
{
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << "CGAL_CSLV relocate" << std::endl;
  std::cout << "HeCv: " << he->curve() << std::endl;
  std::cout << "HeDi: " << he->direction() << std::endl;
#endif

  // We use a constant index-map to prvent the introduction of new entries.
  // When not using a const index-map, erroneous entries might be added!!!
  const Halfedge_indices_map& const_he_indices_table = m_he_indices_table;

  // Go along the boundary of the new face.
  Face_handle new_face = he->face();
  Halfedge_handle curr_he = he;
  const Halfedge_handle invalid_he;

#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << "m_sc_counter: " << m_sc_counter << std::endl;
  std::cout << "m_sc_he_table: " << m_sc_he_table.size() << std::endl;
#endif
  do {
    // We are interested only in halfedges directed from right to left.
    if (curr_he->direction() == ARR_LEFT_TO_RIGHT) {
      curr_he = curr_he->next();
      continue;
    }

    // Get the indices list associated with the current halfedges, representing
    // the halfedges and isolated vertices that "see" it from above.
    const Indices_list& indices_list = const_he_indices_table[curr_he];
    typename Indices_list::const_iterator itr;
    for (itr = indices_list.begin(); itr != indices_list.end(); ++itr) {
      CGAL_assertion(*itr != 0);
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
      std::cout << "itr: " << *itr << std::endl;
#endif

      // In case the current subcurve index does not match a valid entry in
      // m_sc_he_table, we know that this subcurve matches a halfedge that is
      // not yet mapped. This can happen only if this halfedge is he itself.
      // As we know that he lies on the outer CCB of the new face, it is
      // definitely not a hole in the face, therefore we can ignore it.
      if ((*itr > m_sc_counter) || (*itr >= m_sc_he_table.size())) continue;

      Halfedge_handle he_on_face = m_sc_he_table[*itr];
      if (he_on_face == invalid_he) {
        // If the halfedge handle is invalid, then we have an index for an
        // isolated vertex. Move this vertex to the new face, if necessary.
        Vertex_handle v = m_iso_verts_map[*itr];
        CGAL_assertion(v != m_invalid_vertex);
        if (v->face() != new_face) {
          m_arr_access.move_isolated_vertex(v->face(), new_face, v);
        }
      }
      else {
        // If necessary, move the hole that the halfedge belongs to into the
        // new face.
        if ((he_on_face->twin()->face() != new_face) &&
            he_on_face->twin()->is_on_inner_ccb())
        {
          m_arr_access.move_inner_ccb(he_on_face->twin()->face(),
                                      new_face,
                                      he_on_face->twin()->ccb());

          // Perform the relocation process recursively: Namely all holes
          // and isolated vertices that "see" he_on_face from above should also
          // be located inside the new face.
          relocate_in_new_face(he_on_face->twin());
        }
      }
    }
    curr_he = curr_he->next();
  } while(curr_he != he);
}

//-----------------------------------------------------------------------------
// Map the given subcurve index to the given halfedge handle.
//
template <typename Hlpr, typename Vis>
void Arr_construction_ss_visitor<Hlpr, Vis>::
_map_new_halfedge(unsigned int i, Halfedge_handle he)
{
#if CGAL_ARR_CONSTRUCTION_SS_VISITOR_VERBOSE
  std::cout << "map " << i << " to " << he->curve() << " "
            << he->direction() << std::endl;
#endif
  CGAL_assertion(i != 0);
  // Resize the index table if needed.
  if (i >= m_sc_he_table.size()) m_sc_he_table.resize(i+1);

  // Map the index to the given halfedge handle.
  m_sc_he_table[i] = he;
}

} // namespace CGAL

#endif
