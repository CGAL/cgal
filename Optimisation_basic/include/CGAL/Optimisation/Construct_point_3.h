// Copyright (c) 1997-2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_OPTIMISATION_CONSTRUCT_POINT_3_H
#define CGAL_OPTIMISATION_CONSTRUCT_POINT_3_H

#include <CGAL/Point_3.h>
#include <vector>
#include <functional>
#include <iterator>

namespace CGAL {

// Class declaration
// =================
template < class K >
class _Construct_point_3;

// Class interface
// ===============
template < class K_ >
class _Construct_point_3 {
  public:
    // self
    typedef  K_                         K;
    typedef  _Construct_point_3<K>       Self;

    // types
    typedef  typename K::Point_3        Point;

    // creation
    _Construct_point_3( ) { }

    // operations
    template < class InputIterator >
    Point
    operator() ( int, InputIterator first, InputIterator last) const
    {
        InputIterator i(first);
	typename K::RT x = *(i++);
	typename K::RT y = *(i++);
	typename K::RT z = *(i++);
	typedef typename K::Construct_point_3 Construct_point_3;
	Construct_point_3 construct_point_3 = K().construct_point_3_object();
	if (i==last) {
	    return construct_point_3(x,y,z);
	} else {
	    typename K::RT h = *(i++);
	    return construct_point_3(x,y,z,h); 
	}
    }

};

} //namespace CGAL

#endif // CGAL_OPTIMISATION_CONSTRUCT_POINT_3_H

// ===== EOF ==================================================================
