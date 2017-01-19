// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein        <wein@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_DEFAULT_DCEL_H
#define CGAL_ARR_DEFAULT_DCEL_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * The definition of the Arr_default_dcel<Traits> class.
 */

#include <CGAL/Arr_dcel_base.h>

namespace CGAL {

/*! \class
 * The default arrangement DCEL class.
 * The Traits parameters corresponds to a geometric traits class, which 
 * defines the Point_2 and X_monotone_curve_2 types.
 */
template <class Traits_>
class Arr_default_dcel :
  public Arr_dcel_base<Arr_vertex_base<typename Traits_::Point_2>,
		       Arr_halfedge_base<typename Traits_::X_monotone_curve_2>,
		       Arr_face_base>
{
public:

  /*! \struct
   * An auxiliary structure for rebinding the DCEL with a new traits class.
   */
  template<typename T>
  struct rebind
  {
    typedef Arr_default_dcel<T> other;
  };

  /*! Default constructor. */
  Arr_default_dcel()
  {}

  /*! Destructor. */
  virtual ~Arr_default_dcel()
  {}
};

} //namespace CGAL

#endif 
