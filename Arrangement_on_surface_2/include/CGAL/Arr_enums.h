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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s): Ron Wein          <wein@post.tau.ac.il>
//            Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARR_ENUM_H
#define CGAL_ARR_ENUM_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 * Definition of the enumeration types for the arrangement package.
 */

#include <CGAL/basic.h>
#include <CGAL/enum.h>

namespace CGAL {

/*! \enum
 * Selection of a curve end.
 */
enum Arr_curve_end
{
  ARR_MIN_END,
  ARR_MAX_END
};

//! \brief prints curve end (for debugging)
template< class OutputStream >
inline
OutputStream& operator<<(
        OutputStream& os,
        const Arr_curve_end& ce) {
    
    switch(ce) {
    case CGAL::ARR_MIN_END:
        os << "ARR_MIN_END";
        break;
    case CGAL::ARR_MAX_END:
        os << "ARR_MAX_END";
        break;
    default:
        CGAL_error_msg("bogus curve end");
    }
    return os;
}


/*! \enum
 * Indicator whether a halfedge is directed from left to right (from the
 * xy-lexicographically smaller vertex to the larger one), or from right to
 * left.
 */
enum Arr_halfedge_direction
{
  ARR_LEFT_TO_RIGHT = -1,
  ARR_RIGHT_TO_LEFT = 1
};

//! \brief prints halfedge direction (for debugging)
template< class OutputStream >
inline
OutputStream& operator<<(
        OutputStream& os,
        const Arr_halfedge_direction& dir) {
    
    switch(dir) {
    case CGAL::ARR_LEFT_TO_RIGHT:
        os << "ARR_LEFT_TO_RIGHT";
        break;
    case CGAL::ARR_RIGHT_TO_LEFT:
        os << "ARR_RIGHT_TO_LEFT";
        break;
    default:
        CGAL_error_msg("bogus halfedge direction");
    }
    return os;
}

/*! \enum The various surface boundary types.
 * For example:
 * - The plain has unbounded boundaries.
 * - The sphere has 2 contraction points and one identification curve.
 */
enum Arr_boundary_type {
  ARR_OBLIVIOUS = 0,
  ARR_OPEN,
  ARR_CLOSED,
  ARR_CONTRACTION,
  ARR_IDENTIFICATION
};

//! \brief prints boundary type (for debugging)
template< class OutputStream >
inline
OutputStream& operator<<(
        OutputStream& os,
        const Arr_boundary_type& bt) {
    
    switch(bt) {
    case CGAL::ARR_OPEN:
        os << "ARR_OPEN";
        break;
    case CGAL::ARR_CLOSED:
        os << "ARR_CLOSED";
        break;
    case CGAL::ARR_CONTRACTION:
        os << "ARR_CONTRACTION";
        break;
    case CGAL::ARR_IDENTIFICATION:
        os << "ARR_IDENTIFICATION";
        break;
    case CGAL::ARR_OBLIVIOUS:
        os << "ARR_OBLIVIOUS";
        break;
    default:
        CGAL_error_msg("bogus boundary type");
    }
    return os;
}


/*! \enum The various surface parameter space options categorizing the
 * surface range according to the parameter domain.
 */

typedef Box_parameter_space_2 Arr_parameter_space;

const Arr_parameter_space ARR_LEFT_BOUNDARY = LEFT_BOUNDARY;

const Arr_parameter_space ARR_RIGHT_BOUNDARY = RIGHT_BOUNDARY;

const Arr_parameter_space ARR_BOTTOM_BOUNDARY = BOTTOM_BOUNDARY;

const Arr_parameter_space ARR_TOP_BOUNDARY = TOP_BOUNDARY;

const Arr_parameter_space ARR_INTERIOR = INTERIOR;


//! \brief prints parameter space (for debugging)
template< class OutputStream >
inline
OutputStream& operator<<(
        OutputStream& os,
        const Arr_parameter_space& ps) {
  
  switch (::CGAL::get_mode(os)) {
  case ::CGAL::IO::PRETTY:
    switch(ps) {
    case CGAL::ARR_LEFT_BOUNDARY:
      os << "ARR_LEFT_BOUNDARY";
      break;
    case CGAL::ARR_RIGHT_BOUNDARY:
      os << "ARR_RIGHT_BOUNDARY";
      break;
    case CGAL::ARR_BOTTOM_BOUNDARY:
      os << "ARR_BOTTOM_BOUNDARY";
      break;
    case CGAL::ARR_TOP_BOUNDARY:
      os << "ARR_TOP_BOUNDARY";
      break;
    case CGAL::ARR_INTERIOR:
      os << "ARR_INTERIOR";
      break;
    default:
      CGAL_error_msg("bogus parameter space");
    }
    break;
  case ::CGAL::IO::BINARY:
    std::cerr << "BINARY format not yet implemented" << std::endl;
    break;
  default:
    os << static_cast< int >(ps);
  }
  
  return os;
}


//! \brief reads parameter space
template< class InputStream >
inline
InputStream& operator>>(
    InputStream& is,
    Arr_parameter_space& ps) {

  CGAL_precondition(CGAL::is_ascii(is));

  int i;
  is >> i;
  ps = static_cast< Arr_parameter_space >(i);
  
  return is;
  
}


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
