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
// Author(s): Ron Wein          <wein@post.tau.ac.il>
//            Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARRANGEMENT_WITH_HISTORY_2_H
#define CGAL_ARRANGEMENT_WITH_HISTORY_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The header file for the `Arrangement_with_history_2<GeomTraits, Dcel>` class.
 */

#include <utility>

#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arrangement_on_surface_with_history_2.h>
#include <CGAL/Arrangement_2/Arr_default_planar_topology.h>

namespace CGAL {

/*! \class Arrangement_with_history_2
 * The arrangement with history class, representing planar subdivisions
 * induced by a set of arbitrary planar curves and storing the curve history.
 * The `GeomTraits` parameter corresponds to a geometry-traits class that
 * defines the `Point_2`, `X_monotone_curve_2`, and `Curve_2` types and implements the
 * geometric predicates and constructions for the family of curves it defines.
 * The `Dcel` parameter should be a model of the `AosDcelWithRebind` concept and support
 * the basic topological operations on a doubly-connected edge-list.
 */
template <typename GeomTraits_, typename Dcel_ = Arr_default_dcel<GeomTraits_>>
class Arrangement_with_history_2 :
  public Arrangement_on_surface_with_history_2<GeomTraits_,
                                               typename Default_planar_topology<GeomTraits_, Dcel_>::Traits> {
private:
  using Default_topology = Default_planar_topology<GeomTraits_, Dcel_>;
  using Base = Arrangement_on_surface_with_history_2<GeomTraits_, typename Default_topology::Traits>;

public:
  using Geometry_traits_2 = GeomTraits_;
  using Dcel = Dcel_;

  using Point_2 = typename Base::Point_2;
  using X_monotone_curve_2 = typename Base::X_monotone_curve_2;
  using Curve_2 = typename Base::Curve_2;

  using Topology_traits = typename Base::Topology_traits;

  // Type definitions.
  using Vertex = typename Base::Vertex;
  using Halfedge = typename Base::Halfedge;
  using Face = typename Base::Face;
  using Size = typename Base::Size;

  using Vertex_iterator = typename Base::Vertex_iterator;
  using Vertex_const_iterator = typename Base::Vertex_const_iterator;

  using Halfedge_iterator = typename Base::Halfedge_iterator;
  using Halfedge_const_iterator = typename Base::Halfedge_const_iterator;

  using Edge_iterator = typename Base::Edge_iterator;
  using Edge_const_iterator = typename Base::Edge_const_iterator;

  using Face_iterator = typename Base::Face_iterator;
  using Face_const_iterator = typename Base::Face_const_iterator;

  using Halfedge_around_vertex_circulator = typename Base::Halfedge_around_vertex_circulator;
  using Halfedge_around_vertex_const_circulator = typename Base::Halfedge_around_vertex_const_circulator;

  using Ccb_halfedge_circulator = typename Base::Ccb_halfedge_circulator;
  using Ccb_halfedge_const_circulator = typename Base::Ccb_halfedge_const_circulator;

  using Outer_ccb_iterator = typename Base::Outer_ccb_iterator;
  using Outer_ccb_const_iterator = typename Base::Outer_ccb_const_iterator;

  using Inner_ccb_iterator = typename Base::Inner_ccb_iterator;
  using Inner_ccb_const_iterator = typename Base::Inner_ccb_const_iterator;

  using Isolated_vertex_iterator = typename Base::Isolated_vertex_iterator;
  using Isolated_vertex_const_iterator = typename Base::Isolated_vertex_const_iterator;

  using Vertex_handle = typename Base::Vertex_handle;
  using Vertex_const_handle = typename Base::Vertex_const_handle;

  using Halfedge_handle = typename Base::Halfedge_handle;
  using Halfedge_const_handle = typename Base::Halfedge_const_handle;

  using Face_handle = typename Base::Face_handle;
  using Face_const_handle = typename Base::Face_const_handle;

  using Curve_iterator = typename Base::Curve_iterator;
  using Curve_const_iterator = typename Base::Curve_const_iterator;
  using Curve_handle = typename Base::Curve_handle;
  using Curve_const_handle = typename Base::Curve_const_handle;
  using Originating_curve_iterator = typename Base::Originating_curve_iterator;
  using Induced_edge_iterator = typename Base::Induced_edge_iterator;

  /*! a shared pointer to the (immutable) geometry traits. */
  using Shared_geometry_traits = typename Base::Shared_geometry_traits;

  // These types are defined for backward compatibility:
  using Traits_2 = Geometry_traits_2;
  using Hole_iterator = typename Base::Inner_ccb_iterator;
  using Hole_const_iterator = typename Base::Inner_ccb_const_iterator;

private:
  using Self = Arrangement_with_history_2<Geometry_traits_2, Dcel>;

  friend class Arr_accessor<Self>;

public:
  /// \name Constructors.
  //@{

  /*! constructs default. */
  Arrangement_with_history_2() = default;

  /*! constructs copy (from a base arrangement). */
  Arrangement_with_history_2(const Base& base) : Base(base) {}

  /*! constructs given a shared traits object. The arrangement (co-)owns the traits. */
  explicit Arrangement_with_history_2(Shared_geometry_traits tr) : Base(std::move(tr)) {}

  /*! constructs from a traits object (owned by the caller). */
  Arrangement_with_history_2(const Traits_2* tr) : Base(tr) {}
  //@}

  /// \name Assignment functions.
  //@{

  /*! assigns (from a base arrangement). */
  Self& operator=(const Base& base) {
    Base::assign(base);
    return *this;
  }

  /*! assigns an arrangement. */
  void assign(const Base& base) { Base::assign(base); }
  //@}

  /// \name Specialized access methods.
  //@{

  /*! obtains the geometry-traits class (for backward compatibility). */
  const Traits_2* traits() const { return this->geometry_traits(); }

  /*! obtains the number of vertices at infinity. */
  Size number_of_vertices_at_infinity() const {
    // The vertices at infinity are valid, but not concrete:
    return this->topology_traits()->number_of_valid_vertices() - this->topology_traits()->number_of_concrete_vertices();
  }

  /*! obtains the unbounded face (non-const version). */
  Face_handle unbounded_face() { return this->non_const_handle(std::as_const(*this).unbounded_face()); }

  /*! obtains the unbounded face (const version). */
  Face_const_handle unbounded_face() const {
    // The fictitious un_face contains all other valid faces in a single
    // hole inside it. We return a handle to one of its neighboring faces,
    // which is necessarily unbounded.
    const typename Base::DFace* un_face = this->topology_traits()->initial_face();

    if (! un_face->is_fictitious()) return Face_const_handle(un_face);

    const typename Base::DHalfedge* p_he = *(un_face->inner_ccbs_begin());
    const typename Base::DHalfedge* p_opp = p_he->opposite();
    const typename Base::DOuter_ccb* p_oc = p_opp->outer_ccb();
    return Face_const_handle(p_oc->face());
  }

  //@}
};

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
