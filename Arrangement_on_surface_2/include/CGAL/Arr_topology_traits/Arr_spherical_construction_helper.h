// Copyright (c) 2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_SPHERICAL_CONSTRUCTION_HELPER_H
#define CGAL_ARR_SPHERICAL_CONSTRUCTION_HELPER_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Definition of the Arr_spherical_construction_helper class-template.
 */

#include <CGAL/Arr_accessor.h>
#include <CGAL/Unique_hash_map.h>

namespace CGAL {

/*! \class Arr_spherical_construction_helper
 *
 * A helper class for the construction sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <typename GeometryTraits_2, typename Arrangement_, typename Event_,
          typename Subcurve_>
class Arr_spherical_construction_helper {
public:
  typedef GeometryTraits_2                              Geometry_traits_2;
  typedef Arrangement_                                  Arrangement_2;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;
  typedef typename Subcurve::Allocator                  Allocator;

protected:
  typedef Geometry_traits_2                             Gt2;

public:
  typedef typename Gt2::X_monotone_curve_2              X_monotone_curve_2;
  typedef typename Gt2::Point_2                         Point_2;

  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;

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
    typedef Arr_spherical_construction_helper<OtherGeometryTraits_2,
                                              OtherArrangement,
                                              OtherEvent, OtherSubcurve>
                                                        other;
  };

  // The following should be private. It is declared 'protected' as a
  // workaround to a problem with VC. (At least VC 14 exhibits this problem).
  // When declared private, VC claims that Gt2 is private (within
  // Arr_spherical_construction_helper); thus, it cannot be access by
  // Arr_spherical_construction_helper.
protected:
  typedef typename Arrangement_2::Topology_traits       Topology_traits;

  typedef typename Topology_traits::Vertex              DVertex;
  typedef typename Topology_traits::Halfedge            DHalfedge;

  // Data members:

  //! The topology-traits class
  Topology_traits* m_top_traits;

  //! An arrangement accessor
  Arr_accessor<Arrangement_2> m_arr_access;

  //! The unbounded arrangement face
  Face_handle m_spherical_face;

  //! Indices of the curves that "see" the north pole.
  Indices_list m_subcurves_at_nf;

  //! A pointer to a map of halfedges to indices lists
  // (stored in the visitor class)
  Halfedge_indices_map* m_he_ind_map_p;

public:
  /*! Constructor. */
  Arr_spherical_construction_helper(Arrangement_2* arr) :
    m_top_traits(arr->topology_traits()),
    m_arr_access(*arr),
    m_he_ind_map_p(nullptr)
  {}

  /*! Destructor. */
  virtual ~Arr_spherical_construction_helper() {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  virtual void before_sweep()
  { m_spherical_face = Face_handle(m_top_traits->spherical_face()); }

  /*! A notification invoked before the sweep-line starts handling the given
   * event.
   */
  virtual void before_handle_event(Event* event);

  /*! A notification invoked when a new subcurve is created. */
  virtual void add_subcurve(Halfedge_handle , Subcurve* ) { return; }

  /*! Collect a subcurve index that does not see any status-line from below.
   */
  void add_subcurve_in_top_face(unsigned int index)
  { m_subcurves_at_nf.push_back(index); }

  /*! Obtain the indices of the halfedges that "see" the north. */
  Indices_list& halfedge_indices_list() { return m_subcurves_at_nf; }

  /*! A notification invoked before the given event it deallocated. */
  void before_deallocate_event(Event* /* event */) {}
  //@}

  /*! Set the map that maps each halfedge to the list of subcurve indices
   * that "see" the halfedge from below.
   */
  void set_halfedge_indices_map(Halfedge_indices_map& table)
  { m_he_ind_map_p = &table; }

  /*! Determine if we should swap the order of predecessor halfedges when
   * calling insert_at_vertices_ex() .
   */
  bool swap_predecessors(Event* event) const
  {
    // If we insert an edge whose right end lies on the north pole, we have
    // to flip the order of predecessor halfegdes.
    return (event->parameter_space_in_x() == ARR_INTERIOR &&
            event->parameter_space_in_y() == ARR_TOP_BOUNDARY);
  }

  /*! Get the current top face. */
  virtual Face_handle top_face() const { return m_spherical_face; }

  /*! Splice the indices list of the given halfedge, inserting the
   * indices of the halfedges that "see" the north pole.
   */
  void splice_indices_list(Halfedge_handle he)
  {
    CGAL_assertion(m_he_ind_map_p != nullptr);
    Indices_list& list_ref = (*m_he_ind_map_p)[he];
    list_ref.splice(list_ref.end(), m_subcurves_at_nf);
  }
};

/*! A notification invoked before the sweep-line starts handling the given
 * event.
 */
template <typename Traits_, typename Arr_, typename Event_, typename Subcurve_>
void Arr_spherical_construction_helper<Traits_, Arr_, Event_, Subcurve_>::
before_handle_event(Event* event)
{
  // Act according to the boundary type:
  Arr_parameter_space ps_x = event->parameter_space_in_x();
  Arr_parameter_space ps_y = event->parameter_space_in_y();
  if ((ps_x == ARR_INTERIOR) && (ps_y == ARR_INTERIOR)) return;

  if (ps_y == ARR_BOTTOM_BOUNDARY) {
    // Process bootom contraction boundary:
    // The event has only one right curve, as there is exactly one curve
    // incident to an event with boundary conditions.
    CGAL_assertion((event->number_of_left_curves() == 0) &&
                   (event->number_of_right_curves() == 1));
    const X_monotone_curve_2& xc =
      (*(event->right_curves_begin()))->last_curve();

    // If a vertex on the south pole does not exists, create one.
    DVertex* dv = m_top_traits->south_pole();
    Vertex_handle v = (dv) ? Vertex_handle(dv) :
      m_arr_access.create_boundary_vertex(xc, ARR_MIN_END, ps_x, ps_y);
    event->set_vertex_handle(v);
    return;
  }

  if (ps_y == ARR_TOP_BOUNDARY) {
    // Process top contraction boundary:
    // The event has only one left curve, as there is exactly one curve
    // incident to an event with boundary conditions.
    CGAL_assertion((event->number_of_left_curves() == 1) &&
                   (event->number_of_right_curves() == 0));
    const X_monotone_curve_2& xc =
      (*(event->left_curves_begin()))->last_curve();

    DVertex* dv = m_top_traits->north_pole();
    if (dv) {
      event->set_vertex_handle(Vertex_handle(dv));
      DHalfedge* dprev =
        m_top_traits->locate_around_boundary_vertex(m_top_traits->north_pole(),
                                                    xc, ARR_MAX_END,
                                                    ps_x, ps_y);

      if (!dprev) return;
      Halfedge_handle prev = Halfedge_handle(dprev);
      event->set_halfedge_handle(prev);

      // Associate the indices of all subcurves that "see" the top boundary
      // (from below) with the halfedge that (clockwise) succeeds the
      // predecessor of the curve to be inserted around the north pole. The
      // next halfedge is actually the last halfedge that starts at the
      // north pole encountered. Observe that the direction of every
      // halfedge that starts at the north pole is right to left.
      splice_indices_list(prev->next());
      return;
    }

    // We do not have a vertex that corresponds to the north pole.
    // Create one.
    Vertex_handle v =
      m_arr_access.create_boundary_vertex(xc, ARR_MAX_END, ps_x, ps_y);
    event->set_vertex_handle(v);

    // Since this is the first event corresponding to the north pole,
    // the list m_subcurves_at_nf contains all subcurves whose minimal
    // endpoint lies between the curve of discontinuity and the current
    // curve incident to the north pole. In case these subcurves represent
    // holes, these holes should stay in the "north face" that contains the
    // line of discontinuity, and we should not keep track of them in order
    // to later move them to another face.
    m_subcurves_at_nf.clear();
    return;
  }

  if (ps_x == ARR_LEFT_BOUNDARY) {
    // The event has only right curves, as there is exactly one curve
    // incident to an event with boundary conditions.
    CGAL_assertion((event->number_of_left_curves() == 0) &&
                   (event->number_of_right_curves() >= 1));
    const X_monotone_curve_2& xc =
      (*(event->right_curves_begin()))->last_curve();

    // If a vertex on the line of discontinuity does not exists. create one.
    DVertex* dv = m_top_traits->discontinuity_vertex(xc, ARR_MIN_END);
    Vertex_handle v = (dv) ? Vertex_handle(dv) :
      m_arr_access.create_boundary_vertex(xc, ARR_MIN_END, ps_x, ps_y);
    event->set_vertex_handle(v);
    return;
  }

  if (ps_x == ARR_RIGHT_BOUNDARY) {
    // The event has only left curves, as there is exactly one curve
    // incident to an event with boundary conditions.
    CGAL_assertion((event->number_of_left_curves() >= 1) &&
                   (event->number_of_right_curves() == 0));
    const X_monotone_curve_2& xc =
      (*(event->left_curves_begin()))->last_curve();

    // If a vertex on the line of discontinuity does not exists. create one.
    DVertex* dv = m_top_traits->discontinuity_vertex(xc, ARR_MAX_END);
    Vertex_handle v = (dv) ? Vertex_handle(dv) :
      m_arr_access.create_boundary_vertex(xc, ARR_MAX_END, ps_x, ps_y);
    event->set_vertex_handle(v);
    return;
  }
}

} // namespace CGAL

#endif
