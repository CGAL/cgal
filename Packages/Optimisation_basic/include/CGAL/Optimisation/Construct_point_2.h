// Copyright (c) 1997-2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_OPTIMISATION_CONSTRUCT_POINT_2_H
#define CGAL_OPTIMISATION_CONSTRUCT_POINT_2_H

// includes
#ifndef CGAL_POINT_2_H
#  include <CGAL/Point_2.h>
#endif
#ifndef CGAL_PROTECT_VECTOR
#  include <vector>
#  define CGAL_PROTECT_VECTOR
#endif
#ifndef CGAL_PROTECT_FUNCTIONAL
#  include <functional>
#  define CGAL_PROTECT_FUNCTIONAL
#endif
#ifndef CGAL_PROTECT_ITERATOR
#  include <iterator>
#  define CGAL_PROTECT_ITERATOR
#endif

CGAL_BEGIN_NAMESPACE

// Class declaration
// =================
template < class R >
class Construct_point_2;

// Class interface
// ===============
template < class R_ >
class Construct_point_2 {
  public:
    // self
    typedef  R_                         R;
    typedef  Construct_point_2<R>       Self;

    // types
    typedef  typename R::Point_2        Point;

    // creation
    Construct_point_2( ) { }

    // operations
    template < class InputIterator >
    Point
    operator() ( int, InputIterator first, InputIterator last) const
    {
	std::vector<typename R::RT>  coords;
	std::copy( first, last, std::back_inserter( coords));
	if ( coords.size() < 3) {
	    return Point( coords[ 0], coords[ 1]);
	} else {
	    return Point( coords[ 0], coords[ 1], coords[ 2]);
	}
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_OPTIMISATION_CONSTRUCT_POINT_2_H

// ===== EOF ==================================================================
