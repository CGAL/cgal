// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://golubevs@scm.gforge.inria.fr/svn/cgal/branches/experimental-packages/Arrangement_on_surface_2/include/CGAL/Arr_topology_traits/Arr_bounded_planar_construction_helper.h $
// $Id: Arr_bounded_planar_construction_helper.h 34244 2006-09-14 13:14:21Z wein $
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Ron Wein <wein@post.tau.ac.il>

#ifndef CGAL_ARR_BOUNDED_PLANAR_CONSTRUCTION_HELPER_H
#define CGAL_ARR_BOUNDED_PLANAR_CONSTRUCTION_HELPER_H

/*!
 * Definition of the Arr_bounded_planar_construction_helper class-template.
 */

#include <CGAL/Sweep_line_empty_visitor.h>
#include <CGAL/Unique_hash_map.h>

CGAL_BEGIN_NAMESPACE

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

  typedef Sweep_line_empty_visitor<Traits_2,
                                   Subcurve,
                                   Event>              Base_visitor;
public:
 
  /*! Constructor. */
  Arr_bounded_planar_construction_helper(Arrangement_2 *arr=0)
  {}

  /*! Destructor. */
  virtual ~Arr_bounded_planar_construction_helper()
  {}

};

CGAL_END_NAMESPACE

#endif
