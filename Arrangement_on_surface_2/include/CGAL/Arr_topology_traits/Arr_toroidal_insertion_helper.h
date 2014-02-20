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

#ifndef CGAL_ARR_TOROIDAL_INSERTION_HELPER_H
#define CGAL_ARR_TOROIDAL_INSERTION_HELPER_H

/*! \file
 * Definition of the Arr_toroidal_insertion_helper class-template.
 */

#include <CGAL/Sweep_line_2/Arr_construction_sl_visitor.h>
#include <CGAL/Arr_topology_traits/Arr_toroidal_construction_helper.h>

namespace CGAL {

/*! \class Arr_toroidal_insertion_helper
 * A helper class for the insertion sweep-line visitors, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <typename Traits_, typename Arrangement_, typename Event_,
          typename Subcurve_>
class Arr_toroidal_insertion_helper :
  public Arr_toroidal_construction_helper<Traits_, Arrangement_,
                                         Event_, Subcurve_>
{
public:
  typedef Traits_                                       Traits_2;
  typedef Arrangement_                                  Arrangement_2;
  typedef Event_                                        Event;
  typedef Subcurve_                                     Subcurve;

  typedef typename Traits_2::X_monotone_curve_2         X_monotone_curve_2;
  typedef typename Traits_2::Point_2                    Point_2;

  typedef Arr_toroidal_construction_helper<Traits_2, Arrangement_2, Event,
                                            Subcurve>   Base;

  typedef Sweep_line_empty_visitor<Traits_2, Subcurve, Event>
                                                        Base_visitor;

  typedef Arr_toroidal_insertion_helper<Traits_2, Arrangement_2, Event,
                                         Subcurve>      Self;

  typedef Arr_construction_sl_visitor<Self>             Parent_visitor;

  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;

protected:
  typedef typename Arrangement_2::Topology_traits       Topology_traits;
  typedef typename Topology_traits::Vertex              DVertex;

public:
  /*! Constructor */
  Arr_toroidal_insertion_helper(Arrangement_2 *arr) : Base(arr) {}

  /*! Destructor. */
  virtual ~Arr_toroidal_insertion_helper() {}

  /// \name Notification functions.
  //@{

  /* A notification issued before the sweep process starts. */
  virtual void before_sweep() { return; }

  /*! A notification invoked before the sweep-line starts handling a given
   * event.
   */
  virtual void before_handle_event(Event* event) { return; }

  /*! A notification invoked when a new subcurve is created. */
  virtual void add_subcurve(Halfedge_handle he, Subcurve* sc) { return; }

  //!@}

}; // Arr_toroidal_insertion_helper

} //namespace CGAL

#endif // CGAL_ARR_TOROIDAL_INSERTION_HELPER_H
