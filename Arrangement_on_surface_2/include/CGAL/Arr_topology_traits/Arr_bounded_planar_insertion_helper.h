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

#ifndef CGAL_ARR_BOUNDED_PLANAR_INSERTION_HELPER_H
#define CGAL_ARR_BOUNDED_PLANAR_INSERTION_HELPER_H

/*!
 * Definition of the Arr_bounded_planar_insertion_helper class-template.
 */

#include <CGAL/Sweep_line_2/Arr_construction_sl_visitor.h>
#include <CGAL/Arr_topology_traits/Arr_bounded_planar_construction_helper.h>

namespace CGAL {

/*! \class Arr_bounded_planar_insertion_helper
 * A helper class for the insertion sweep-line visitors, suitable
 * for an Arrangement_on_surface_2 instantiated with a topology-traits class
 * for bounded curves in the plane.
 */
template <class Traits_, class Arrangement_, class Event_, class Subcurve_> 
class Arr_bounded_planar_insertion_helper :
  public Arr_bounded_planar_construction_helper<Traits_, Arrangement_,
                                                Event_, Subcurve_>
{
public:

  typedef Traits_                                      Traits_2;
  typedef Arrangement_                                 Arrangement_2;
  typedef Event_                                       Event;
  typedef Subcurve_                                    Subcurve;

  typedef typename Traits_2::X_monotone_curve_2        X_monotone_curve_2;
  typedef typename Traits_2::Point_2                   Point_2;

  typedef Arr_bounded_planar_construction_helper<Traits_2,
                                             Arrangement_2,
                                             Event,
                                             Subcurve> Base;

  typedef Sweep_line_empty_visitor<Traits_2,
                                   Subcurve,
                                   Event>              Base_visitor;

  typedef Arr_bounded_planar_insertion_helper<Traits_2,
                                          Arrangement_2,
                                          Event,
                                          Subcurve>    Self;

  typedef Arr_construction_sl_visitor<Self>            Parent_visitor;

  typedef typename Arrangement_2::Face_handle          Face_handle;

  typedef typename Base::Indices_list                  Indices_list;
  typedef typename Base::Halfedge_indices_map          Halfedge_indices_map;

public:
 
  /*! Constructor. */
  Arr_bounded_planar_insertion_helper (Arrangement_2 *arr) :
    Base (arr)
  {}
};

} //namespace CGAL

#endif
