// Copyright (c) 2007 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
//

#ifndef CGAL_ARR_SPHERICAL_VERT_DECOMP_HELPER_H
#define CGAL_ARR_SPHERICAL_VERT_DECOMP_HELPER_H

/*! \file
 * Definition of the Arr_spherical_vert_decomp_helper class-template.
 */

CGAL_BEGIN_NAMESPACE

#include <CGAL/Sweep_line_empty_visitor.h>

/*! \class Arr_spherical_vert_decomp_helper
 * A helper class for the vertical decomposition sweep-line visitor, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <class Traits_, class Arrangement_>
class Arr_spherical_vert_decomp_helper
{
public:

  typedef Traits_                                      Traits_2;
  typedef Arrangement_                                 Arrangement_2;

  typedef typename Arrangement_2::Face_const_handle    Face_const_handle;

  typedef Sweep_line_empty_visitor<Traits_2>           Base_visitor;
  typedef typename Base_visitor::Event                 Event;
  typedef typename Base_visitor::Subcurve              Subcurve;

protected:

  typedef typename Arrangement_2::Topology_traits      Topology_traits;

  const Topology_traits  *m_top_traits;     // The topology traits.
  Face_const_handle       m_north_face;     // Current north face.

public:

  /*! Constructor.
   * \param arr The arrangement.
   */
  Arr_spherical_vert_decomp_helper(const Arrangement_2 *arr) :
    m_top_traits(arr->topology_traits())
  {}

  /// \name Notification functions.
  //@{

  /*! A notification issued before the sweep process starts. */
  void before_sweep()
  {
    // Get the spherical face of the arrangement.
    m_north_face = Face_const_handle(m_top_traits->spherical_face());
    return;
  }

  /*!
   * A notification invoked after the sweep-line finishes handling the given
   * event.
   */
  void after_handle_event(Event * event)
  {
    return;
  }
  //@}

  /*! Get the current top object. */
  CGAL::Object top_object () const
  {
    // Wrap the north face with a CGAL object.
    return (CGAL::make_object (m_north_face));
  }

  /*! Get the current bottom object. */
  CGAL::Object bottom_object () const
  {
    // Wrap the north face with a CGAL object.
    return (CGAL::make_object (m_north_face));
  }
};

CGAL_END_NAMESPACE

#endif
