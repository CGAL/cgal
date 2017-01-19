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

#ifndef CGAL_ARR_BOUNDED_PLANAR_CONSTRUCTION_HELPER_H
#define CGAL_ARR_BOUNDED_PLANAR_CONSTRUCTION_HELPER_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*!
 * Definition of the Arr_bounded_planar_construction_helper class-template.
 */

#include <CGAL/Sweep_line_empty_visitor.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

/*! \class Arr_bounded_planar_construction_helper
 * A helper class for the construction sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <class Traits_, class Arrangement_, class Event_, class Subcurve_> 
class Arr_bounded_planar_construction_helper
{
public:

  typedef Traits_                                      Traits_2;
  typedef Arrangement_                                 Arrangement_2;
  typedef Event_                                       Event;
  typedef Subcurve_                                    Subcurve;

  typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename Traits_2::Point_2                   Point_2;

  typedef Sweep_line_empty_visitor<Traits_2,
                                   Subcurve,
                                   Event>              Base_visitor;

  typedef typename Arrangement_2::Face_handle          Face_handle;
  typedef typename Arrangement_2::Halfedge_handle      Halfedge_handle;
  
  typedef typename Subcurve::Halfedge_indices_list     Indices_list;
  typedef Unique_hash_map<Halfedge_handle, 
                          Indices_list>                Halfedge_indices_map;

protected:

  typedef typename Arrangement_2::Topology_traits      Topology_traits;

  // Data members:
  Topology_traits         *m_top_traits;  // The topology-traits class.
  Face_handle              m_unb_face;    // The unbounded arrangement face.

  Indices_list m_emptylist;

public:
 
  /*! Constructor. */
  Arr_bounded_planar_construction_helper (Arrangement_2 *arr) :
    m_top_traits (arr->topology_traits())
  {}

  /*! Destructor. */
  virtual ~Arr_bounded_planar_construction_helper()
  {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  virtual void before_sweep ()
  {
    // Get the unbounded face.
    m_unb_face = Face_handle (m_top_traits->unbounded_face());
  }

  /*!
   * A notification invoked before the sweep-line starts handling the given
   * event.
   */
  virtual void before_handle_event (Event* /* event */)
  {
    return;
  }

  /*! A notification invoked when a new subcurve is created. */
  virtual void add_subcurve (Halfedge_handle /* he */, Subcurve* /* sc */)
  {
    return;
  }

  Indices_list& halfedge_indices_list() {
    return m_emptylist;
  }

  /*! Collect a subcurve index that does not see any status-line from below. */
  void add_subcurve_in_top_face (unsigned int /* index */)
  {
    return;
  }

  /*! A notification invoked before the given event it deallocated. */
  void before_deallocate_event (Event* /* event */)
  {
    return;
  }
  //@} 
  
  /*!
   * Set the map that maps each halfedge to the list of subcurve indices
   * that "see" the halfedge from below.
   */
  void set_halfedge_indices_map (Halfedge_indices_map& /* table */)
  {
    return;
  }

  /*!
   * Determine if we should swap the order of predecessor halfedges when
   * calling insert_at_vertices_ex() .
   */
  bool swap_predecessors (Event* /* event */) const
  {
    // In the bounded case the order of the predecessor is always correct
    // and there is no need to swap them.
    return (false);
  }

  /*! Get the current top face. */
  Face_handle top_face () const
  {
    return (m_unb_face);
  }
};

} //namespace CGAL

#endif
