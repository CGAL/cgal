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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_BOUNDED_PLANAR_CONSTRUCTION_HELPER_H
#define CGAL_ARR_BOUNDED_PLANAR_CONSTRUCTION_HELPER_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_bounded_planar_construction_helper class-template.
 */

#include <CGAL/Unique_hash_map.h>

namespace CGAL {

/*! \class Arr_bounded_planar_construction_helper
 *
 * A helper class for the construction sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <typename GeometryTraits_2, typename Arrangement_, typename Event_,
          typename Subcurve_>
class Arr_bounded_planar_construction_helper {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Arrangement_                                  Arrangement_2;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;
  typedef typename Subcurve::Allocator                  Allocator;

  // The following should be private. It is declared 'protected' as a
  // workaround to a problem with VC. (At least VC 14 exhibits this problem).
  // When declared private, VC claims that Gt2 is private (within
  // Arr_bounded_planar_construction_helper); thus, it cannot be access by
  // Arr_bounded_planar_construction_helper.
protected:
  typedef Geometry_traits_2                             Gt2;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;

  typedef typename Subcurve::Halfedge_indices_list      Indices_list;
  typedef Unique_hash_map<Halfedge_handle, Indices_list>
                                                        Halfedge_indices_map;

  /*! \struct rebind
   * An auxiliary structure for rebinding the helper with a new types.
   * Mainly used to rebind the geometry-traits type and a new type that derives
   * from the old one.
   */
  template <typename OtherGeometryTraits_2, typename OtherArrangement,
            typename OtherEvent, typename OtherSubcurve>
  struct rebind {
    typedef Arr_bounded_planar_construction_helper<OtherGeometryTraits_2,
                                                   OtherArrangement,
                                                   OtherEvent, OtherSubcurve>
                                                        other;
  };

protected:
  typedef typename Arrangement_2::Topology_traits       Topology_traits;

  // Data members:
  Topology_traits* m_top_traits;        // The topology-traits class.
  Face_handle m_unb_face;               // The unbounded arrangement face.

  Indices_list m_emptylist;

public:
  /*! Constructor. */
  Arr_bounded_planar_construction_helper(Arrangement_2* arr) :
    m_top_traits(arr->topology_traits())
  {}

  /*! Destructor. */
  virtual ~Arr_bounded_planar_construction_helper()
  {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  virtual void before_sweep()
  {
    // Get the unbounded face.
    m_unb_face = Face_handle(m_top_traits->unbounded_face());
  }

  /*! A notification invoked before the sweep-line starts handling the given
   * event.
   */
  virtual void before_handle_event(Event* /* event */) { return; }

  /*! A notification invoked when a new subcurve is created. */
  virtual void add_subcurve(Halfedge_handle /* he */, Subcurve* /* sc */)
  { return; }

  Indices_list& halfedge_indices_list() { return m_emptylist; }

  /*! Collect a subcurve index that does not see any status-line from below. */
  void add_subcurve_in_top_face(unsigned int /* index */) { return; }

  /*! A notification invoked before the given event it deallocated. */
  void before_deallocate_event(Event* /* event */) { return; }
  //@}

  /*! Set the map that maps each halfedge to the list of subcurve indices
   * that "see" the halfedge from below.
   */
  void set_halfedge_indices_map(Halfedge_indices_map& /* table */)
  { return; }

  /*! Determine if we should swap the order of predecessor halfedges when
   * calling insert_at_vertices_ex() .
   */
  bool swap_predecessors(Event* /* event */) const
  {
    // In the bounded case the order of the predecessor is always correct
    // and there is no need to swap them.
    return false;
  }

  /*! Get the current top face. */
  Face_handle top_face() const { return m_unb_face; }
};

} // namespace CGAL

#endif
