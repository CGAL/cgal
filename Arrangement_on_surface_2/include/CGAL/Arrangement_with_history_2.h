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

#ifndef CGAL_ARRANGEMENT_WITH_HISTORY_2_H
#define CGAL_ARRANGEMENT_WITH_HISTORY_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * The header file for the Arrangement_with_history_2<Traits,Dcel> class.
 */

#include <CGAL/Arrangement_on_surface_with_history_2.h>
#include <CGAL/Arrangement_2/Arr_default_planar_topology.h>

namespace CGAL {

/*! \class Arrangement_on_history_2
 * The arrangement with history class, representing planar subdivisions
 * induced by a set of arbitrary planar curves and storing the curve history.
 * The GeomTraits parameter corresponds to a geometry-traits class that
 * defines the Point_2 and X_monotone_curve_2 types and implements the
 * geometric predicates and constructions for the family of curves it defines.
 * The Dcel parameter should be a model of the ArrDcel concept and support
 * the basic topological operations on a doubly-connected edge-list.
 */
template <class GeomTraits_,
          class Dcel_ = Arr_default_dcel<GeomTraits_> >
class Arrangement_with_history_2 :
  public Arrangement_on_surface_with_history_2
    <GeomTraits_,
     typename Default_planar_topology<GeomTraits_, Dcel_>::Traits>
{
private:

  typedef Default_planar_topology<GeomTraits_, Dcel_>     Default_topology;
  typedef Arrangement_on_surface_with_history_2<GeomTraits_,
             typename Default_topology::Traits>           Base;

public:
  typedef GeomTraits_                                     Geometry_traits_2;
  typedef Dcel_                                           Dcel;

  typedef typename Base::Point_2                          Point_2;
  typedef typename Base::X_monotone_curve_2               X_monotone_curve_2;
  typedef typename Base::Curve_2                          Curve_2;

  typedef typename Base::Topology_traits                  Topology_traits;

  // Type definitions.
  typedef typename Base::Vertex                   Vertex;
  typedef typename Base::Halfedge                 Halfedge;
  typedef typename Base::Face                     Face;
  typedef typename Base::Size                     Size;

  typedef typename Base::Vertex_iterator          Vertex_iterator;
  typedef typename Base::Vertex_const_iterator    Vertex_const_iterator;

  typedef typename Base::Halfedge_iterator        Halfedge_iterator;
  typedef typename Base::Halfedge_const_iterator  Halfedge_const_iterator;

  typedef typename Base::Edge_iterator            Edge_iterator;
  typedef typename Base::Edge_const_iterator      Edge_const_iterator;

  typedef typename Base::Face_iterator            Face_iterator;
  typedef typename Base::Face_const_iterator      Face_const_iterator;

  typedef typename Base::Halfedge_around_vertex_circulator
                                     Halfedge_around_vertex_circulator;
  typedef typename Base::Halfedge_around_vertex_const_circulator
                                     Halfedge_around_vertex_const_circulator;

  typedef typename Base::Ccb_halfedge_circulator   Ccb_halfedge_circulator;
  typedef typename Base::Ccb_halfedge_const_circulator
                                             Ccb_halfedge_const_circulator;

  typedef typename Base::Outer_ccb_iterator        Outer_ccb_iterator;
  typedef typename Base::Outer_ccb_const_iterator  Outer_ccb_const_iterator;

  typedef typename Base::Inner_ccb_iterator        Inner_ccb_iterator;
  typedef typename Base::Inner_ccb_const_iterator  Inner_ccb_const_iterator;

  typedef typename Base::Isolated_vertex_iterator  Isolated_vertex_iterator;
  typedef typename Base::Isolated_vertex_const_iterator
                                             Isolated_vertex_const_iterator;

  typedef typename Base::Vertex_handle             Vertex_handle;
  typedef typename Base::Vertex_const_handle       Vertex_const_handle;

  typedef typename Base::Halfedge_handle           Halfedge_handle;
  typedef typename Base::Halfedge_const_handle     Halfedge_const_handle;

  typedef typename Base::Face_handle               Face_handle;
  typedef typename Base::Face_const_handle         Face_const_handle;

  typedef typename Base::Curve_iterator             Curve_iterator;
  typedef typename Base::Curve_const_iterator       Curve_const_iterator;
  typedef typename Base::Curve_handle               Curve_handle;
  typedef typename Base::Curve_const_handle         Curve_const_handle;
  typedef typename Base::Originating_curve_iterator Originating_curve_iterator;
  typedef typename Base::Induced_edge_iterator      Induced_edge_iterator;

  // These types are defined for backward compatibility:
  typedef Geometry_traits_2                        Traits_2;
  typedef typename Base::Inner_ccb_iterator        Hole_iterator;
  typedef typename Base::Inner_ccb_const_iterator  Hole_const_iterator;

private:
  typedef Arrangement_with_history_2<Geometry_traits_2, Dcel>   Self;

  friend class Arr_observer<Self>;
  friend class Arr_accessor<Self>;

public:

  /// \name Constructors.
  //@{

  /*! Default constructor. */
  Arrangement_with_history_2 () :
    Base ()
  {}

  /*! Copy constructor (from a base arrangement). */
  Arrangement_with_history_2 (const Base & base) :
    Base (base)
  {}

  /*! Constructor given a traits object. */
  Arrangement_with_history_2 (const Traits_2 * tr) :
    Base (tr)
  {}
  //@}

  /// \name Assignment functions.
  //@{

  /*! Assignment operator (from a base arrangement). */
  Self& operator= (const Base & base)
  {
    Base::assign (base);
    return (*this);
  }

  /*! Assign an arrangement. */
  void assign (const Base & base)
  {
    Base::assign (base);
    return;
  }
  //@}

  ///! \name Specialized access methods.
  //@{

  /*! Get the geometry-traits class (for backward compatibility). */
  const Traits_2 * traits () const
  {
    return (this->geometry_traits());
  }

  /*! Get the number of vertices at infinity. */
  Size number_of_vertices_at_infinity () const
  {
    // The vertices at infinity are valid, but not concrete:
    return (this->topology_traits()->number_of_valid_vertices() -
            this->topology_traits()->number_of_concrete_vertices());
  }

  /*! Get the number of unbounded faces. */
  Size number_of_unbounded_faces () const
  {
    typename Base::Face_const_iterator  fit = this->faces_begin();
    typename Base::Face_const_iterator  end = this->faces_end();
    Size                                n_unb = 0;

    while (fit != end)
    {
      if (fit->is_unbounded())
        n_unb++;
      ++fit;
    }

    return (n_unb);
  }

  /*! Get the unbounded face (non-const version). */
  Face_handle unbounded_face ()
  {
    // The fictitious un_face contains all other valid faces in a single
    // hole inside it. We return a handle to one of its neighboring faces,
    // which is necessarily unbounded.
    typename Base::DFace      *un_face =
      const_cast<typename Base::DFace*>(this->topology_traits()->
                                        initial_face());

    if (! un_face->is_fictitious())
      return (Face_handle (un_face));

    typename Base::DHalfedge  *p_he = *(un_face->inner_ccbs_begin());
    typename Base::DHalfedge  *p_opp = p_he->opposite();
    typename Base::DOuter_ccb *p_oc = p_opp->outer_ccb();

    return (Face_handle (p_oc->face()));
  }

  /*! Get the unbounded face (const version). */
  Face_const_handle unbounded_face () const
  {
    // The fictitious un_face contains all other valid faces in a single
    // hole inside it. We return a handle to one of its neighboring faces,
    // which is necessarily unbounded.
    const typename Base::DFace      *un_face =
                                   this->topology_traits()->initial_face();

    if (! un_face->is_fictitious())
      return (Face_const_handle (un_face));

    const typename Base::DHalfedge  *p_he = *(un_face->inner_ccbs_begin());
    const typename Base::DHalfedge  *p_opp = p_he->opposite();
    const typename Base::DOuter_ccb *p_oc = p_opp->outer_ccb();

    return (Face_const_handle (p_oc->face()));
  }
  //@}

protected:

  /// \name Managing and notifying the arrangement observers.
  //@{
  typedef Arr_observer<Self>                      Observer;

  /*!
   * Register a new observer (so it starts receiving notifications).
   * \param p_obs A pointer to the observer object.
   */
  void _register_observer (Observer *p_obs)
  {
    Base::_register_observer ((typename Base::Observer*)p_obs);
    return;
  }

  /*!
   * Unregister a new observer (so it stops receiving notifications).
   * \param p_obs A pointer to the observer object.
   * \return Whether the observer was successfully unregistered.
   */
  bool _unregister_observer (Observer *p_obs)
  {
    return (Base::_unregister_observer ((typename Base::Observer*)p_obs));
  }
  //@}

};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
