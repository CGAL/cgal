// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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

#ifndef CGAL_ARR_CONSTRUCTION_SL_VISITOR_H
#define CGAL_ARR_CONSTRUCTION_SL_VISITOR_H

#ifndef CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE
#define CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE 0
#endif

#ifndef CGAL_NEW_FACE_SPLIT_STRATEGY 
#define CGAL_NEW_FACE_SPLIT_STRATEGY 0
#endif

/*!
 * Definition of the Arr_construction_sl_visitor class-template.
 */

#include <CGAL/Arr_accessor.h>
#include <CGAL/Unique_hash_map.h>
#include <vector>

namespace CGAL {

/*! \struct Integer_hash_function
 * An auxiliary hash functor for integers.
 */
struct Integer_hash_function 
{
  typedef std::size_t result_type;

  std::size_t operator() (unsigned int i) const 
  { 
    return i;
  }  
};

/*! \class Arr_construction_sl_visitor
 * A sweep-line visitor for constructing an arrangement embedded on a surface.
 */
template <class Helper_> 
class Arr_construction_sl_visitor :
  public Helper_::Base_visitor
{
public:

  typedef Helper_                                      Helper;

  typedef typename Helper::Traits_2                    Traits_2;
  typedef typename Helper::Arrangement_2               Arrangement_2;
  typedef typename Helper::Base_visitor                Base;
  typedef typename Helper::Event                       Event;
  typedef typename Helper::Subcurve                    Subcurve;

  typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename Traits_2::Point_2                   Point_2;

protected:

  typedef typename Arrangement_2::Topology_traits      Topology_traits;
  typedef typename Arrangement_2::Vertex_handle        Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle      Halfedge_handle;
  typedef typename Arrangement_2::Face_handle          Face_handle;

  typedef typename Base::Event_subcurve_iterator       Event_subcurve_iterator;
  typedef typename Base::Event_subcurve_reverse_iterator
                                               Event_subcurve_reverse_iterator;
  typedef typename Base::Status_line_iterator          Status_line_iterator;

  typedef typename Helper::Indices_list                Indices_list;
  typedef typename Helper::Halfedge_indices_map        Halfedge_indices_map;
  typedef Unique_hash_map<unsigned int,
                          Vertex_handle,
                          Integer_hash_function>       Iso_vertices_map;

protected:

  Helper                   m_helper;      // The helper class.

  Arrangement_2           *m_arr;         // The arrangement we construct.
  Topology_traits         *m_top_traits;  // The topology-traits class.
  Arr_accessor<Arrangement_2>
                           m_arr_access;  // An arrangement accessor.

  unsigned int             m_sc_counter;  // Counter for subcurves that may
                                          // represent a hole (the upper
                                          // subcurves that emarge from event
                                          // points with only right curves). 

  std::vector<Halfedge_handle>
                           m_sc_he_table; // A table that maps a subcurve
                                          // index to its halfedge handle,
                                          // directed from right to left.

  Iso_vertices_map         m_iso_verts_map; // Maps an index to the isolated
                                            // vertex.

  Halfedge_indices_map     m_he_indices_table; // Maps each halfdge to the
                                               // indices of subcurves that
                                               // lies below it.

  const Vertex_handle      m_invalid_vertex;   // An invalid vertex handle.

public:
 
  /*! Constructor. */
  Arr_construction_sl_visitor (Arrangement_2 *arr) :
    m_helper (arr),
    m_arr (arr),
    m_top_traits (arr->topology_traits()),
    m_arr_access (*arr),
    m_sc_counter (0),
    m_sc_he_table (1),
    m_invalid_vertex ()
  {
    m_helper.set_halfedge_indices_map (m_he_indices_table);
  }

  /*! Destructor. */
  virtual ~Arr_construction_sl_visitor()
  {}

  /// \name Sweep-line notifications.
  //@{

  /* A notification issued before the sweep process starts. */
  inline void before_sweep ();

  /*!
   * A notification invoked before the sweep-line starts handling the given
   * event.
   */
  inline void before_handle_event (Event* event);

  /*!
   * A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  bool after_handle_event (Event* event, Status_line_iterator iter, bool flag);

  /*! A notification invoked when a new subcurve is created. */
  void add_subcurve (const X_monotone_curve_2& cv, Subcurve* sc);
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
  insert_in_face_interior (const X_monotone_curve_2& cv, Subcurve* sc);

  /*!
   * Insert the given subcurve given its left end-vertex.
   * \param cv The geometric entity.
   * \param prev The predecessor halfedge around the left vertex.
   * \param sc The sweep-line subcurve information.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle
  insert_from_left_vertex (const X_monotone_curve_2& cv, Halfedge_handle he,
                           Subcurve* sc);

  /*!
   * Insert the given subcurve given its right end-vertex.
   * \param cv The geometric entity.
   * \param prev The predecessor halfedge around the right vertex.
   * \param sc The sweep-line subcurve information.
   * \return A handle to the inserted halfedge.
   */
  virtual Halfedge_handle
  insert_from_right_vertex (const X_monotone_curve_2& cv,
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
  virtual Halfedge_handle insert_at_vertices (const X_monotone_curve_2& cv,
                                              Halfedge_handle prev1,
                                              Halfedge_handle prev2,
                                              Subcurve* sc,
                                              bool &new_face_created);

  /*!
   * Insert an isolated vertex into the arrangement.
   * \param pt The point associated with the vertex.
   * \param iter The location of the corresponding event in the status line.
   * \return A handle to the inserted vertex.
   */
  virtual Vertex_handle insert_isolated_vertex (const Point_2& pt,
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
  Event* last_event_on_subcurve (Subcurve* sc)
  {
    return (reinterpret_cast<Event*>((sc)->last_event()));
  }

private:

  /// \name Auxiliary functions.
  //@{

  /*!
   * Cast a Traits::Point_2 object into an Arrangement_2::Point_2 object.
   * These two types may not be the same when the addition visitor inherits
   * from this base class.
   */
  inline const typename Arrangement_2::Point_2& _point (const Point_2& p) const
  {
    return (static_cast<const typename Arrangement_2::Point_2&> (p));
  }

  /*!
   * Cast a Traits::X_monotone_curve_2 object into an
   * Arrangement_2::X_monotone_curve_2 object.
   * These two types may not be the same when the addition visitor inherits
   * from this base class.
   */
  inline const typename Arrangement_2::X_monotone_curve_2&
  _curve (const X_monotone_curve_2& cv) const
  {
    return (static_cast<const typename Arrangement_2::X_monotone_curve_2&>(cv));
  }

  /*! Map the given subcurve index to the given halfedge handle. */
  void _map_new_halfedge (unsigned int i, Halfedge_handle he);
  //@}
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
template <class Hlpr> 
void Arr_construction_sl_visitor<Hlpr>::before_sweep ()
{
  // We just have to notify the helper that the sweep process now starts.
  m_helper.before_sweep();

  return;
}

//-----------------------------------------------------------------------------
// A notification invoked before the sweep-line starts handling the given
// event.
//
template <class Hlpr> 
void Arr_construction_sl_visitor<Hlpr>::before_handle_event (Event* event)
{
  // We just have to notify the helper class on the event.
  m_helper.before_handle_event (event);
#if 0
  std::cout << "CGAL_CSLV before_handle_event" << std::endl;
#endif
  return;
}

//-----------------------------------------------------------------------------
// A notification invoked after the sweep-line finishes handling the given
// event.
//
template <class Hlpr> 
bool Arr_construction_sl_visitor<Hlpr>::after_handle_event
    (Event* event, Status_line_iterator iter, bool /* flag */)
{
#if 0
    std::cout << "CGAL_CSLV after_handle_event" << std::endl;
#endif

  // Check if the event represents an isolated vertex.
  if (!event->has_left_curves() && !event->has_right_curves())
  {
    // There are no incident subcurves, so this event is an isolated vertex.
    // We map the current index to this vertex, and add this index to the
    // indices list of the curve the vertex "sees" from below.
    Vertex_handle  v = insert_isolated_vertex(event->point(), iter);

    m_sc_counter++;
    m_iso_verts_map[m_sc_counter] = v;
    _map_new_halfedge(m_sc_counter, Halfedge_handle());

    if (iter != this->status_line_end())
    {
      // The isolated vertex "sees" the subcurve of the given position from
      // below.
      Subcurve *sc_above = *iter;
      sc_above->add_halfedge_index(m_sc_counter);
    }
    else
    {
      // The vertex is not located below any valid curve, so we use the helper
      // class to mark that this index should belong to the current top face.
#if 0
        std::cout << "CGAL_CSLV adding a " << m_sc_counter << std::endl;
#endif
      m_helper.add_subcurve_in_top_face (m_sc_counter);
    }

    // The event can now be deallocated.
    return (true);
  }

  // Check if the event has only incident subcurves from its right.
  if (!event->has_left_curves() && !event->is_on_boundary())
  {
    CGAL_assertion(event->has_right_curves());

    // In case of a finite event that has no incident left curves, it is
    // associated with a point that may be the leftmost one in a hole.
    // We give index to the topmost subcurve from the right, and add this
    // vertex indices list of the curve the event "sees" from below.
    m_sc_counter++;
    (*(event->right_curves_rbegin()))->set_index (m_sc_counter);

    if (iter != this->status_line_end())
    {
      // The vertex "sees" the subcurve of the given position from below.
      Subcurve *sc_above = *iter;
      sc_above->add_halfedge_index(m_sc_counter);
    }
    else
    {
      // The vertex is not located below any valid curve, so we use the helper
      // class to mark that this index should belong to the current top face.
#if 0
        std::cout << "CGAL_CSLV adding b " << m_sc_counter << std::endl;
#endif
      m_helper.add_subcurve_in_top_face (m_sc_counter);
    }
  }
  
  // Set the last event of all left subcurve (thus, this event corresponds
  // to their right endpoint).
  Event_subcurve_iterator  left_it;
  for(left_it = event->left_curves_begin();
      left_it != event->left_curves_end();
      ++left_it)
  {
    (*left_it)->set_last_event(event);
  }
  
  // In case there are no right subcurves, the event can be deallocated.
  if(event->number_of_right_curves() == 0)
  {
    // Inform the helper class that the event will soon be deallocated.
    m_helper.before_deallocate_event (event);
    return (true);
  }

  // Mark that all right subcurves incident to the current event are not
  // in the arrangement yet.
  event->init_subcurve_in_arrangement_flags (event->number_of_right_curves());
  
  // Set the last event of all right subcurve (thus, this event corresponds
  // to their left endpoint).
  Event_subcurve_iterator  right_it;
  for(right_it = event->right_curves_begin();
      right_it != event->right_curves_end();
      ++right_it)
  {
    (*right_it)->set_last_event(event);
  }

  // Mark that the event cannot be deallocated just yet.
  return (false); 
}

//-----------------------------------------------------------------------------
// A notification invoked when a new subcurve is created.
//
template <class Hlpr> 
void Arr_construction_sl_visitor<Hlpr>::
add_subcurve (const X_monotone_curve_2& cv, Subcurve* sc)
{
#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE
    std::cout << "CGAL_CSLV add_subcurve: " << cv << std::endl;
#endif

  // Obtain all information to perform the insertion of the subcurve into
  // the arrangement.
  Event           *last_event = last_event_on_subcurve(sc);
  Halfedge_handle  res; 
  Halfedge_handle  he_right = this->current_event()->halfedge_handle();
  Halfedge_handle  he_left = last_event->halfedge_handle();
  const int jump = last_event->compute_halfedge_jump_count(sc);

  const Halfedge_handle  invalid_he;

#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE 
  if (last_event->is_closed()) {
      std::cout << "CGAL_CSLG lastevent: " << last_event->point() << std::endl;
  }
  if (he_left != invalid_he) {
      if (!he_left->is_fictitious()) {
          std::cout << "he_leftcv  : " << he_left->curve() << std::endl;
      } else {
          std::cout << "he_left    : fictitious" << std::endl;
      }
      std::cout << "he_leftdir : " << he_left->direction() << std::endl;
      std::cout << "he_leftfac : " << &(*he_left->face()) << std::endl;
  } else {
      std::cout << "he_left    : invalid" << std::endl;
  }
  if (he_right != invalid_he) {
      if (!he_right->is_fictitious()) {
          std::cout << "he_rightcv : " << he_right->curve() << std::endl;
      } else {
          std::cout << "he_right   : fictitious" << std::endl;
      }
      std::cout << "he_rightdir: " << he_right->direction() << std::endl;
      std::cout << "he_rightfac: " << &(*he_right->face()) << std::endl;
  } else {
      std::cout << "he_right   : invalid" << std::endl;
  }
#endif

  // Check if the previous event on the curve is not in the arrangement yet.
  if (he_left == invalid_he) 
  {
    // We do not have a handle from the previous insert.
    if (he_right != invalid_he)
    {
      // We have a handle from the current event, representing the right end
      // of the subcurve - use it to insert the subcurve.
#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE 
        std::cout << "CGAL_CSLV call insert_from_right_vertex" << std::endl;
#endif
      res = this->insert_from_right_vertex(cv, he_right, sc);
    }
    else
    {
      // We do not have handles for any of the curve end, so we insert it in
      // the interior of a face.
#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE 
        std::cout << "CGAL_CSLV call insert_in_face_interior" << std::endl;
#endif
      res = this->insert_in_face_interior(cv, sc);
    }
  } 
  else 
  {
    // The previous event on the curve is already in the arrangement,
    // therefore  we use it to insert the subcurve.
    // First, we skip some halfedges around the left vertex to get the true
    // predecessor halfedge for the insertion.
    for (int i = 0; i < jump; i++) {
      he_left = (he_left->next())->twin();
    }
#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE
      if (jump != 0) {
          std::cout << "CGAL_CSLV JUMP: " << jump << std::endl;
          if (!he_left->is_fictitious()) {
              std::cout << "he_leftcv  : " << he_left->curve() << std::endl;
          } else {
              std::cout << "he_left    : fictitious" << std::endl;
          }
          std::cout << "he_leftdir : " << he_left->direction() << std::endl;
          std::cout << "he_leftfac : " << &(*he_left->face()) << std::endl;
      }
#endif
      
    if (he_right != invalid_he) 
    {
      CGAL_assertion (he_left->face() == he_right->face());
      
      // We also have a handle for the current event, representing the right
      // vertex of the subcurve. We insert the subcurve using the two
      // predecessor halfedges.       
      bool   dummy;
      
#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE 
      std::cout << "CGAL_CSLV call insert_at_vertices" << std::endl;
#endif
      res = this->insert_at_vertices (cv, he_right, he_left, sc, dummy);
    }
    else
    {
      // We only have a handle for the predecessor halfedge of the left end
      // of the subcurve - use it to insert the subcurve.
#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE 
        std::cout << "CGAL_CSLV call insert_from_left_vertex" << std::endl;
#endif
      res = this->insert_from_left_vertex (cv, he_left, sc);
    }
  }

  // Make sure that res is a halfedge that is always directed from left to
  // right (thus its twin is directed from right to left).
  if (res->direction() != ARR_LEFT_TO_RIGHT)
    res = res->twin();

  // Update the last event with the inserted halfegde (if necessary)
  // and check if we have to update the auxiliary information on the location 
  // of holes.
  if (last_event->number_of_left_curves() == 0 &&  
      last_event->is_curve_largest((Subcurve*)sc))
  {
      if (last_event->vertex_handle() == m_invalid_vertex)
          last_event->set_halfedge_handle(res->twin());
    
    // If sc has valid index, insert its index to m_sc_he_table.
    if(sc->has_valid_index())
    {
      CGAL_assertion(res->twin()->direction() == ARR_RIGHT_TO_LEFT);
      _map_new_halfedge (sc->index(), res->twin());
    }
  }

  // Update the halfedge handle associated with the current event.
  if (this->current_event()->vertex_handle() == m_invalid_vertex)
      this->current_event()->set_halfedge_handle(res);

  // In case the event has no more right subcurves associated with it, we can
  // deallocate it. Note that we inform the helper class before deallocating
  // the event.
  if (last_event->dec_right_curves_counter() == 0)
  {
    m_helper.before_deallocate_event (last_event);
    this->deallocate_event (last_event);
  }

  // Clear the list of indices of the subcurve.
  sc->clear_halfedge_indices();
  
  return;
}

//-----------------------------------------------------------------------------
// Insert the given subcurve in the interior of an arrangement face.
//
template <class Hlpr>
typename Arr_construction_sl_visitor<Hlpr>::Halfedge_handle
Arr_construction_sl_visitor<Hlpr>::
insert_in_face_interior (const X_monotone_curve_2& cv, Subcurve* sc)
{
#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE
    std::cout << "CGAL_CSLV insert_in_face_interior\ncurve: " << cv << std::endl;
#endif

  // Check if the vertex to be associated with the left end of the curve has
  // already been created.
  Event         *last_event = last_event_on_subcurve(sc);
  Vertex_handle  v1 = last_event->vertex_handle();
  bool           create_v1 = false;

  if (v1 == m_invalid_vertex)
  {
    // Mark that we should create the vertex v1 later on (if we created it
    // now, and ended up calling insert_from_right_vertex(), this vertex
    // would be constructed twice!)
    create_v1 = true;
  }
  else if (v1->degree() > 0)
  {
    // In this case the left vertex v1 is a boundary vertex which already has
    // some incident halfedges. We look for the predecessor halfedge and
    // and insert the curve from this left vertex.
    Arr_parameter_space   bx = last_event->parameter_space_in_x();
    Arr_parameter_space   by = last_event->parameter_space_in_y();

    CGAL_assertion (bx != ARR_INTERIOR || by != ARR_INTERIOR);

    Halfedge_handle l_prev = Halfedge_handle
      (m_top_traits->locate_around_boundary_vertex (&(*v1), _curve(cv),
                                                    ARR_MIN_END, bx, by));
    
    return (this->insert_from_left_vertex (cv, l_prev, sc));
  }

  // Check if the vertex to be associated with the right end of the curve has
  // already been created.
  Event         *curr_event = this->current_event();
  Vertex_handle  v2 = curr_event->vertex_handle();

  if (v2 == m_invalid_vertex)
  {
    // Create the vertex to be associated with the right end of the curve.
    v2 = m_arr_access.create_vertex (_point (curr_event->point()));
  }
  else if (v2->degree() > 0)
  {
    // In this case the right vertex v2 is a boundary vertex which already has
    // some incident halfedges. We look for the predecessor halfedge and
    // and insert the curve from this right vertex.
    Arr_parameter_space   bx = curr_event->parameter_space_in_x();
    Arr_parameter_space   by = curr_event->parameter_space_in_y();

    CGAL_assertion (bx != ARR_INTERIOR || by != ARR_INTERIOR);

    Halfedge_handle r_prev = Halfedge_handle
      (m_top_traits->locate_around_boundary_vertex (&(*v2), _curve(cv),
                                                    ARR_MAX_END, bx, by));
    
    return (this->insert_from_right_vertex (cv, r_prev, sc));
  }

  // If necessary, create the vertex to be associated with the left end
  // of the curve.
  if (create_v1)
    v1 = m_arr_access.create_vertex (_point (last_event->point()));

  // Perform the insertion between the two (currently isolated) vertices in
  // the interior of the current top face, as given by the helper class.
  Halfedge_handle  res =
    m_arr_access.insert_in_face_interior_ex (_curve(cv), m_helper.top_face(),
                                             v1, v2, SMALLER);

  // Map the new halfedge to the indices list of all subcurves that lie
  // below it.
  if (sc->has_halfedge_indices())
  {
    CGAL_assertion(res->twin()->direction() == ARR_RIGHT_TO_LEFT);
    Indices_list& list_ref = m_he_indices_table[res->twin()];
    list_ref.clear();
    list_ref.splice(list_ref.end(), sc->halfedge_indices_list());
  }

  // Notify the helper on the creation of the new halfedge.
  m_helper.add_subcurve (res, sc);

  return (res);
}

//-----------------------------------------------------------------------------
// Insert the given subcurve using its two end-vertices.
//
template <class Hlpr>
typename Arr_construction_sl_visitor<Hlpr>::Halfedge_handle
Arr_construction_sl_visitor<Hlpr>::insert_at_vertices
    (const X_monotone_curve_2& cv,
     Halfedge_handle prev1,
     Halfedge_handle prev2,
     Subcurve* sc,
     bool& new_face_created)
{

#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE
    std::cout << "CGAL_CSLV insert_at_vertices:\ncurve:" << cv << std::endl;
    if (!prev1->is_fictitious()) {
        std::cout << "prev1cv   : " << prev1->curve() << std::endl;
    } else {
        std::cout << "prev1     : fictitious" << std::endl;
    }
    std::cout << "prev1dir  : " << prev1->direction() << std::endl;
    std::cout << "prev1fac  : " << &(*prev1->face()) << std::endl;
    if (!prev2->is_fictitious()) {
        std::cout << "prev2cv   : " << prev2->curve() << std::endl;
    } else {
        std::cout << "prev2     : fictitious" << std::endl;
    }
    std::cout << "prev2dir  : " << prev2->direction() << std::endl;
    std::cout << "prev2fac  : " << &(*prev2->face()) << std::endl;
#endif

  // Use the helper class to determine whether the order of predecessor
  // halfedges should be swaped, to that the edge directed from prev1->target()
  // to prev2->target() is incident to the new face (in case a new face is
  // created).
  Halfedge_handle   res;
#if CGAL_NEW_FACE_SPLIT_STRATEGY
  // EBEB new strategy for splitting faces. The member also allows
  // to decide which prev_i will be on the new outer CCB (if decision needed)
  // Here: We use it to decide "swapping", which is actually the decision
  //       whether prev1 or prev2 is on the new outer CCB ;-)
  std::pair< bool, bool > update(false, false);
#if 0
  if ((prev1->is_on_inner_ccb() && prev1->is_on_inner_ccb() &&
       prev1->inner_ccb() == prev2->inner_ccb()) 
      ||
      (!prev1->is_on_inner_ccb() && prev1->is_on_inner_ccb())) {
#else
      // TODO improve this code!
  Halfedge_handle curr1 = prev1->next();
  bool found2 = false;
  while (curr1 != prev1) {
      if (curr1 == prev2) {
          found2 = true;
      }
      curr1 = curr1->next();
  }
  Halfedge_handle curr2 = prev2->next();
  bool found1 = false;
  while (curr2 != prev2) {
      if (curr2 == prev1) {
          found1 = true;
      }
      curr2 = curr2->next();
  }
  if (found1 && found2) {
#endif
      update = 
          m_top_traits->face_update_upon_edge_insertion(
                  &(*prev1), &(*prev2), cv
          );
  }
  const bool swap_preds = update.second;
  // TODO propagate update.first to _insert_at_vertices!
#else
  const bool        swap_preds = 
    m_helper.swap_predecessors (this->current_event());
#endif

  res = (! swap_preds) ?
    // usually prev1 is outer of new split face (it it exists)
    m_arr_access.insert_at_vertices_ex (_curve(cv), prev1, prev2,
                                        LARGER, new_face_created) :
    // if swapping prev2 will becomd outer of new split face (it it exists)
    m_arr_access.insert_at_vertices_ex (_curve(cv), prev2, prev1,
                                        SMALLER, new_face_created);
  
  // Map the new halfedge to the indices list of all subcurves that lie
  // below it.
  if (sc->has_halfedge_indices())
  {
    Halfedge_handle he = res;
    
    if (swap_preds)
      he = he->twin();
    CGAL_assertion(he->direction() == ARR_RIGHT_TO_LEFT);

    Indices_list& list_ref = m_he_indices_table[he];
    list_ref.clear();
    list_ref.splice(list_ref.end(), sc->halfedge_indices_list());
  }

  // Notify the helper on the creation of the new halfedge.
  // Note that we do this before trying to relocate holes in the new
  // face (if one is created).
  m_helper.add_subcurve (res, sc);
  
  if (new_face_created)
  {
    // EBEB: Fixed by checking whether at least one of 
    // EBEB: res + res->twin() lies on a inner ccb
    if (res->is_on_inner_ccb() || res->twin()->is_on_inner_ccb()) {
      // In case a new face has been created (pointed by the new halfedge
      // we obtained), we have to examine the holes and isolated vertices
      // in the existing face (pointed by the twin halfedge) and relocate
      // the relevant features in the new face.
      CGAL_assertion(res->face() != res->twin()->face());
      this->relocate_in_new_face (res);
    }
  }

  return (res);
}

//-----------------------------------------------------------------------------
// Insert the given subcurve from a vertex that corresponds to its right end.
//
template <class Hlpr>
typename Arr_construction_sl_visitor<Hlpr>::Halfedge_handle
Arr_construction_sl_visitor<Hlpr>::insert_from_right_vertex
    (const X_monotone_curve_2& cv,
     Halfedge_handle prev,
     Subcurve* sc)
{
#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE
    std::cout << "CGAL_CSLV insert_from_right_vertex:\ncurve:" << cv << std::endl;
    if (!prev->is_fictitious()) {
        std::cout << "prevcv    : " << prev->curve() << std::endl;
    } else {
        std::cout << "prev      : fictitious" << std::endl;
    }
    std::cout << "prevdir   : " << prev->direction() << std::endl;
    std::cout << "prevfac   : " << &(*prev->face()) << std::endl;
#endif

  // Check if the vertex to be associated with the left end of the curve has
  // already been created.
  Event         *last_event = last_event_on_subcurve(sc);
  Vertex_handle  v = last_event->vertex_handle();

  if (v == m_invalid_vertex)
  {
    // Create the vertex to be associated with the left end of the curve.
    v = m_arr_access.create_vertex (_point (last_event->point()));
  }
  else if (v->degree() > 0)
  {
    // In this case the left vertex v is a boundary vertex which already has
    // some incident halfedges. We look for the predecessor halfedge and
    // and insert the curve between two existing vertices.
    Arr_parameter_space   bx = last_event->parameter_space_in_x();
    Arr_parameter_space   by = last_event->parameter_space_in_y();

    CGAL_assertion (bx != ARR_INTERIOR || by != ARR_INTERIOR);

    Halfedge_handle l_prev = Halfedge_handle
      (m_top_traits->locate_around_boundary_vertex (&(*v), _curve(cv),
                                                    ARR_MIN_END, bx, by));
    bool            dummy;

    return (this->insert_at_vertices (cv, prev, l_prev, sc, dummy));
  }

  // Insert the curve given its left vertex and the predecessor around the
  // right vertex.
  Halfedge_handle  res =
    m_arr_access.insert_from_vertex_ex (_curve(cv), prev, v, LARGER);

  // Map the new halfedge to the indices list of all subcurves that lie
  // below it.
  if (sc->has_halfedge_indices())
  {
    CGAL_assertion(res->direction() == ARR_RIGHT_TO_LEFT);
    Indices_list& list_ref = m_he_indices_table[res];
    list_ref.clear();
    list_ref.splice(list_ref.end(), sc->halfedge_indices_list());
  }

  // Notify the helper on the creation of the new halfedge.
  m_helper.add_subcurve (res, sc);

  return (res);
}

//-----------------------------------------------------------------------------
// Insert the given subcurve from a vertex that corresponds to its left end.
//
template <class Hlpr>
typename Arr_construction_sl_visitor<Hlpr>::Halfedge_handle
Arr_construction_sl_visitor<Hlpr>::
insert_from_left_vertex (const X_monotone_curve_2& cv,
                         Halfedge_handle prev,
                         Subcurve* sc)
{
#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE
    std::cout << "CGAL_CSLV insert_from_left_vertex:\ncurve:" << cv << std::endl;
    if (!prev->is_fictitious()) {
        std::cout << "prevcv    : " << prev->curve() << std::endl;
    } else {
        std::cout << "prev      : fictitious" << std::endl;
    }
    std::cout << "prevdir   : " << prev->direction() << std::endl;
    std::cout << "prevfac   : " << &(*prev->face()) << std::endl;
#endif
    
  // Check if the vertex to be associated with the right end of the curve has
  // already been created.
  Event         *curr_event = this->current_event();
  Vertex_handle  v = curr_event->vertex_handle();

  if (v == m_invalid_vertex)
  {
    // Create the vertex to be associated with the right end of the curve.
    v = m_arr_access.create_vertex (_point (curr_event->point()));
  }
  else if (v->degree() > 0)
  {
    // In this case the left vertex v is a boundary vertex which already has
    // some incident halfedges. We look for the predecessor halfedge and
    // and insert the curve from this right vertex.
    Arr_parameter_space   bx = curr_event->parameter_space_in_x();
    Arr_parameter_space   by = curr_event->parameter_space_in_y();

    CGAL_assertion (bx != ARR_INTERIOR || by != ARR_INTERIOR);

    Halfedge_handle r_prev = Halfedge_handle
      (m_top_traits->locate_around_boundary_vertex (&(*v), _curve(cv),
                                                    ARR_MAX_END, bx, by)
        );
    bool            dummy;

    return (this->insert_at_vertices (cv, r_prev, prev, sc, dummy));
  }

  // Insert the curve given its right vertex and the predecessor around the
  // left vertex.
  Halfedge_handle  res =
    m_arr_access.insert_from_vertex_ex (_curve(cv), prev, v, SMALLER);

  // Map the new halfedge to the indices list of all subcurves that lie
  // below it.
  if (sc->has_halfedge_indices())
  {
    CGAL_assertion(res->twin()->direction() == ARR_RIGHT_TO_LEFT);
    Indices_list& list_ref = m_he_indices_table[res->twin()];
    list_ref.clear();
    list_ref.splice(list_ref.end(), sc->halfedge_indices_list());
  }

  // Notify the helper on the creation of the new halfedge.
  m_helper.add_subcurve (res, sc);

  return (res);
}

//-----------------------------------------------------------------------------
// Insert an isolated vertex into the arrangement.
//
template <class Hlpr>
typename Arr_construction_sl_visitor<Hlpr>::Vertex_handle
Arr_construction_sl_visitor<Hlpr>::
insert_isolated_vertex (const Point_2& pt,
                        Status_line_iterator /* iter */)
{
#if CGAL_ARR_CONSTRUCTION_SL_VISITOR_VERBOSE
    std::cout << "CGAL_CSLV insert_isolated_vertex:\npoint:" << pt << std::endl;
#endif

  // Insert the isolated vertex in the interior of the current top face, as
  // given by the helper class.
  return (m_arr->insert_in_face_interior (_point(pt),
                                          m_helper.top_face()));
}

//-----------------------------------------------------------------------------
// Reloacte holes and isolated vertices inside a newly created face.
//
template <class Hlpr>
void Arr_construction_sl_visitor<Hlpr>::
relocate_in_new_face (Halfedge_handle he)
{
#if 0
    std::cout << "CGAL_CSLV relocate" << std::endl;
    std::cout << "HeCv: " << he->curve() << std::endl;
    std::cout << "HeDi: " << he->direction() << std::endl;
#endif

  // We use a constant indices map so no new entries are added there.
  const Halfedge_indices_map& const_he_indices_table = m_he_indices_table;

  // Go along the boundary of the new face.
  Face_handle                 new_face = he->face();
  Halfedge_handle             curr_he = he;
  const Halfedge_handle       invalid_he;
  Vertex_handle               v;

  do
  {
    // We are intreseted only in halfedges directed from right to left.
    if (curr_he->direction() == ARR_LEFT_TO_RIGHT)
    {
      curr_he = curr_he->next();
      continue;
    }
    
    // Get the indices list associated with the current halfedges, representing
    // the halfedges and isolated vertices that "see" it from above.
    const Indices_list& indices_list = const_he_indices_table[curr_he];
    typename Indices_list::const_iterator itr;

    for (itr = indices_list.begin(); itr != indices_list.end(); ++itr)
    {
      CGAL_assertion(*itr != 0);
#if 0
      std::cout << "itr: " << *itr << std::endl;
      std::cout << "m_sc_counter: " << m_sc_counter << std::endl;
      std::cout << "m_sc_he_table: " << m_sc_he_table.size() << std::endl;
#endif

      // In case the current subcurve index does not match a valid entry in 
      // m_sc_he_table, we know that this subcurve matches a halfedge that is
      // not yet mapped. This can happen only if this halfedge is he itself.
      // As we know that he lies on the outer CCB of the new face, it is
      // definately not a hole in the face, therefore we can ignore it.
      if (*itr > m_sc_counter || *itr >= m_sc_he_table.size())
        continue; 
      
      Halfedge_handle he_on_face = m_sc_he_table[*itr];
          
      if(he_on_face == invalid_he)
      {
        // If the halfedge handle is invalis, then we have an index for an
        // isolated vertex. Move this vertex to the new face, if necessary.
        v = m_iso_verts_map[*itr];
        
        CGAL_assertion(v != m_invalid_vertex);
        if (v->face() != new_face)
        {
          m_arr_access.move_isolated_vertex (v->face(),
                                             new_face,
                                             v);
        }
      }
      else
      {
        // If necessary, move the hole that the halfedge belongs to into the
        // new face.
        if (he_on_face->twin()->face() != new_face &&
            he_on_face->twin()->is_on_inner_ccb())
        {
          m_arr_access.move_inner_ccb (he_on_face->twin()->face(),
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

  return;
}

//-----------------------------------------------------------------------------
// Map the given subcurve index to the given halfedge handle.
//
template <class Hlpr>
void Arr_construction_sl_visitor<Hlpr>::
_map_new_halfedge (unsigned int i, Halfedge_handle he)
{
  CGAL_assertion (i != 0);
#if 0
  std::cout << "map " << i << " to " << he->curve() << " " 
            << he->direction() << std::endl;
#endif
  if(i >= m_sc_he_table.size())
    // Resize the index table if we reached it capacity.
    m_sc_he_table.resize(2*i);
  
  // Map the index to the given halfedge handle.
  m_sc_he_table[i] = he;
  return;
}

} //namespace CGAL

#endif
