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

#ifndef CGAL_OPTIMISATION_CONSTRUCT_POINT_D_H
#define CGAL_OPTIMISATION_CONSTRUCT_POINT_D_H

// includes
#  include <CGAL/Kernel_d/Interface_classes.h>
#  include <CGAL/Kernel_d/Point_d.h>

CGAL_BEGIN_NAMESPACE

// Class declaration
// =================
template < class R >
class Construct_point_d;

// Class interface
// ===============
template < class R_ >
class Construct_point_d {
  public:
    // self
    typedef  R_                         R;
    typedef  Construct_point_d<R>       Self;

    // types
    typedef  typename R::Point_d        Point;

    // creation
    Construct_point_d( ) { }

    // operations
    template < class InputIterator >
    Point
    operator() ( int d, InputIterator first, InputIterator last) const
    {
	return Point( d, first, last);
    }
};

CGAL_END_NAMESPACE

#endif // CGAL_OPTIMISATION_CONSTRUCT_POINT_D_H

// ===== EOF ==================================================================
