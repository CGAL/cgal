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
// Author(s)     : Ron Wein <wein@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>
//                 Eric Berberich <eric.berberich@cgal.org>

#ifndef CGAL_ARR_TOROIDAL_VERT_DECOMP_HELPER_H
#define CGAL_ARR_TOROIDAL_VERT_DECOMP_HELPER_H

/*! \file
 * Definition of the Arr_toroidal_vert_decomp_helper class-template.
 */

namespace CGAL {

#include <CGAL/Sweep_line_empty_visitor.h>

/*! \class Arr_toroidal_vert_decomp_helper
 * A helper class for the vertical decomposition sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class.
 */
template <class Traits_, class Arrangement_>
class Arr_toroidal_vert_decomp_helper
{
public:

  typedef Traits_                                      Traits_2;
  typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
  typedef Arrangement_                                 Arrangement_2;

  typedef typename Arrangement_2::Face_const_handle    Face_const_handle;
  typedef typename Arrangement_2::Vertex_const_handle  Vertex_const_handle;

  typedef Sweep_line_empty_visitor<Traits_2>           Base_visitor;
  typedef typename Base_visitor::Event                 Event;
  typedef typename Base_visitor::Subcurve              Subcurve;

protected:

  typedef typename Arrangement_2::Topology_traits      Topology_traits;

  const Topology_traits  *m_top_traits;        // The topology traits.

public:

  /*! Constructor.
   * \param arr The arrangement.
   */
  Arr_toroidal_vert_decomp_helper(const Arrangement_2 *arr) :
    m_top_traits(arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /*! A notification issued before the sweep process starts. */
  void before_sweep();

  /*!
   * A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  void after_handle_event(Event * event);
  //@}

  /*! Get the current top object. */
  CGAL::Object top_object () const;

  /*! Get the current bottom object. */
  CGAL::Object bottom_object () const;
};

//-----------------------------------------------------------------------------
// Memeber-function definitions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// A notification issued before the sweep process starts.
//
template <class Tr, class Arr>
void Arr_toroidal_vert_decomp_helper<Tr, Arr>::before_sweep()
{
  CGAL_error();
  /* dummy implementation */ return;
}

//-----------------------------------------------------------------------------
// A notification invoked after the sweep-line finishes handling the given
// event.
///
template <class Tr, class Arr>
void
Arr_toroidal_vert_decomp_helper<Tr, Arr>::after_handle_event (Event *event)
{
  CGAL_error();
  /* dummy implementation */  return;
}

//-----------------------------------------------------------------------------
// Get the current top object.
///
template <class Tr, class Arr>
CGAL::Object
Arr_toroidal_vert_decomp_helper<Tr, Arr>::top_object () const
{
  CGAL_error();
  /* dummy implementation */ return CGAL::Object();
}

//-----------------------------------------------------------------------------
// Get the current bottom object.
///
template <class Tr, class Arr>
CGAL::Object
Arr_toroidal_vert_decomp_helper<Tr, Arr>::bottom_object () const
{
  CGAL_error();
  /* dummy implementation */ return CGAL::Object();

}

} //namespace CGAL

#endif // CGAL_ARR_TOROIDAL_VERT_DECOMP_HELPER_H
