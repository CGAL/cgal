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

#ifndef CGAL_OPTIMISATION_ACCESS_COORDINATES_BEGIN_D_H
#define CGAL_OPTIMISATION_ACCESS_COORDINATES_BEGIN_D_H


namespace CGAL {

class Cartesian_tag;
class Homogeneous_tag;

// Class declaration
// =================
template < class R_ >
class Access_coordinates_begin_d;

// Class interface
// ===============
template < class R_ >
class Access_coordinates_begin_d {
  public:
    // self
    typedef  R_                         R;
    typedef  Access_coordinates_begin_d<R>
                                        Self;

    // types
    typedef  typename R::Point_d        Point;
    typedef  typename Point::Homogeneous_const_iterator Coordinate_iterator;

    // unary function class types
    typedef  Coordinate_iterator        result_type;
    typedef  Point                      argument_type;

    // creation
    Access_coordinates_begin_d( ) { }

    // operations
private:
    Coordinate_iterator
    access( const Point& p, const Cartesian_tag&) const 
    { return p.cartesian_begin(); }
  
    Coordinate_iterator
    access( const Point& p, const Homogeneous_tag&) const 
    { return p.homogeneous_begin(); }
  
  
public:
    Coordinate_iterator
    operator() ( const Point& p) const { 
      return p.homogeneous_begin();
    }
};

} //namespace CGAL

#endif // CGAL_OPTIMISATION_ACCESS_COORDINATES_BEGIN_D_H

// ===== EOF ==================================================================
