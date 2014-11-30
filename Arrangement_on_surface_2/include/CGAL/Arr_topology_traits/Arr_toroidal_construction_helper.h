// Copyright (c) 2007,2009,2010,2011,2013,2014 Max-Planck-Institute Saarbruecken (Germany), Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//                 Eric Berberich <eric.berberich@cgal.org>

#ifndef CGAL_ARR_TOROIDAL_CONSTRUCTION_HELPER_H
#define CGAL_ARR_TOROIDAL_CONSTRUCTION_HELPER_H

/*! \file
 * Definition of the Arr_toroidal_construction_helper class-template.
 */

#include <CGAL/Sweep_line_empty_visitor.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

/*! \class Arr_toroidal_construction_helper
 * A template helper class for the construction sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class.
 */
template <typename Traits_, typename Arrangement_, typename Event_,
          typename Subcurve_>
class Arr_toroidal_construction_helper {
public:
  typedef Traits_                                         Traits_2;
  typedef Arrangement_                                    Arrangement_2;
  typedef Event_                                          Event;
  typedef Subcurve_                                       Subcurve;

  typedef typename Traits_2::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits_2::Point_2                      Point_2;

  typedef Sweep_line_empty_visitor<Traits_2, Subcurve, Event>
                                                          Base_visitor;

  typedef typename Arrangement_2::Vertex_handle           Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement_2::Face_handle             Face_handle;

  typedef typename Subcurve::Halfedge_indices_list        Indices_list;
  typedef Unique_hash_map<Halfedge_handle, Indices_list>  Halfedge_indices_map;

protected:
  typedef typename Arrangement_2::Topology_traits         Topology_traits;

  typedef typename Topology_traits::Vertex                DVertex;
  typedef typename Topology_traits::Halfedge              DHalfedge;

  // Data members:

  //! The topology-traits class
  Topology_traits* m_top_traits;

  //! An arrangement accessor
  Arr_accessor<Arrangement_2> m_arr_access;

  //! The unbounded arrangement face
  Face_handle m_top_face;

  //! Indices of the curves that "see" the northeast point.
  Indices_list m_subcurves_at_tf;

  //! A pointer to a map of halfedges to indices lists
  // (stored in the visitor class)
  Halfedge_indices_map* m_he_ind_map_p;

public:
  /*! Constructor. */
  Arr_toroidal_construction_helper(Arrangement_2* arr) :
    m_top_traits(arr->topology_traits()),
    m_arr_access(*arr),
    m_he_ind_map_p(NULL)
  {}

  /*! Destructor. */
  virtual ~Arr_toroidal_construction_helper() {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  virtual void before_sweep() {
    m_top_face = Face_handle(m_top_traits->reference_face());
  }

  /*! A notification invoked before the sweep-line starts handling the given
   * event.
   */
  virtual void before_handle_event(Event* event) { return; }

  /*! A notification invoked when a new subcurve is created. */
  virtual void add_subcurve(Halfedge_handle , Subcurve* ) { return; }

  /*! Collect a subcurve index that does not see any status-line from below.
   */
  void add_subcurve_in_top_face(unsigned int index) {
    m_subcurves_at_tf.push_back(index);
  }

  /*! Obtain the indices of the halfedges that "see" the north(east). */
  Indices_list& halfedge_indices_list() { return m_subcurves_at_tf; }

  /*! A notification invoked before the given event it deallocated. */
  void before_deallocate_event(Event* /* event */) { return; }

  //@}

  /*! Set the map that maps each halfedge to the list of subcurve indices
   * that "see" the halfedge from below.
   */
  void set_halfedge_indices_map(Halfedge_indices_map& table) {
    m_he_ind_map_p = &table;
  }

  /*! Determine if we should swap the order of predecessor halfedges when
   * calling insert_at_vertices_ex() .
   */
  bool swap_predecessors(Event* event) const {
    CGAL_error();
    /* dummy implementation */ return false;
  }

  /*! Get the current top face. */
  virtual Face_handle top_face() const {
    return m_top_face;
  }

  /*! Splice the indices list of the given halfedge, inserting the
   * indices of the halfedges that "see" the northeast end.
   */
  void splice_indices_list(Halfedge_handle he)
  {
    CGAL_assertion(m_he_ind_map_p != NULL);
    Indices_list& list_ref = (*m_he_ind_map_p)[he];
    list_ref.splice(list_ref.end(), m_subcurves_at_tf);
  }

}; // Arr_toroidal_construction_helper

} //namespace CGAL

#endif // CGAL_ARR_TOROIDAL_CONSTRUCTION_HELPER_H
