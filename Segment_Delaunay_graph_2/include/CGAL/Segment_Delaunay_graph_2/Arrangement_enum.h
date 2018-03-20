// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARRANGEMENT_ENUM_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARRANGEMENT_ENUM_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>
#include <iostream>


namespace CGAL {

namespace SegmentDelaunayGraph_2 {

namespace Internal {
  struct Arrangement_enum {
    enum Arrangement_type {
      DISJOINT = 0, // obvious
      TOUCH_1, // (p1,p2) and q, and p1 and q are identical
      TOUCH_2, // (p1,p2) and q, and p2 and q are identical
      TOUCH_11, // (p1,p2), (q1,q2), and p1, q1 are identical
      TOUCH_12, // (p1,p2), (q1,q2), and p1, q2 are identical
      TOUCH_21, // (p1,p2), (q1,q2), and p2, q1 are identical
      TOUCH_22, // (p1,p2), (q1,q2), and p2, q2 are identical
      CROSSING, // two segments intersecting at interior points
      IDENTICAL, // either two segments or two points that are identical
      INTERIOR_1, // (p1,p2) and (q1,q2), and q1, q2 are interior
                  // points of (p1,p2)
      INTERIOR_2, // (p1,p2) and (q1,q2), and p1, p2 are interior
                  // points of (q1,q2)
      INTERIOR,  // (p1,p2) and q, and q is an interior point of (p1,p2)
      TOUCH_11_INTERIOR_1, // (p1,p2) and (q1,q2), and p1, q1 are
			   // identical and q2 is an interior point of (p1,p2)
      TOUCH_11_INTERIOR_2, // (p1,p2) and (q1,q2), and p1, q1 are
			   // identical and p2 is an interior point of
			   // (q1,q2)
      TOUCH_12_INTERIOR_1, // (p1,p2) and (q1,q2), and p1, q2 are
                           // identical and q1 is an interior point of (p1,p2)
      TOUCH_12_INTERIOR_2, // (p1,p2) and (q1,q2), and p1, q2 are
			   // identical and p2 is an interior point of (q1,q2)
      TOUCH_21_INTERIOR_1, // (p1,p2) and (q1,q2), and p2, q1 are
			   // identical and q2 is an interior point of (p1,p2)
      TOUCH_21_INTERIOR_2, // (p1,p2) and (q1,q2), and p2, q1 are
  			   // identical and p1 is an interior point of (q1,q2)
      TOUCH_22_INTERIOR_1, // (p1,p2) and (q1,q2), and p2, q2 are
			   // identical and q1 is an interior point of (p1,p2)
      TOUCH_22_INTERIOR_2, // (p1,p2) and (q1,q2), and p2, q2 are
		           // identical and p1 is an interior point of (q1,q2)
      OVERLAPPING_11, // (p1,p2) and (q1,q2), and (p1,q1) is the overlap
      OVERLAPPING_12, // (p1,p2) and (q1,q2), and (p1,q1) is the overlap
      OVERLAPPING_21, // (p1,p2) and (q1,q2), and (p2,q1) is the overlap
      OVERLAPPING_22, // (p1,p2) and (q1,q2), and (p2,q2) is the overlap
      TOUCH_INTERIOR_12, // (p1,p2) and (q1,q2) and p1 is an interior
                         //  point of (q1,q2)
      TOUCH_INTERIOR_22, // (p1,p2) and (q1,q2) and p2 is an interior
                         //  point of (q1,q2)
      TOUCH_INTERIOR_11, // (p1,p2) and (q1,q2) and q1 is an interior
                         //  point of (p1,p2)
      TOUCH_INTERIOR_21  // (p1,p2) and (q1,q2) and q2 is an interior
                         //  point of (p1,p2)
    };


    static Arrangement_type opposite(const Arrangement_type& at) {
      // this returns the result if we swap the order of the arguments...
      if ( at == TOUCH_12 ) {
	return TOUCH_21;
      } else if ( at == TOUCH_21 ) {
	return TOUCH_12;
      } else if ( at == INTERIOR_1 ) {
	return INTERIOR_2;
      } else if ( at == INTERIOR_2 ) {
	return INTERIOR_1;
      } else if ( at == TOUCH_11_INTERIOR_1 ) {
	return TOUCH_11_INTERIOR_2;
      } else if ( at == TOUCH_11_INTERIOR_2 ) {
	return TOUCH_11_INTERIOR_1;
      } else if ( at == TOUCH_12_INTERIOR_1 ) {
	return TOUCH_21_INTERIOR_2;
      } else if ( at == TOUCH_12_INTERIOR_2 ) {
	return TOUCH_21_INTERIOR_1;
      } else if ( at == TOUCH_21_INTERIOR_1 ) {
	return TOUCH_12_INTERIOR_2;
      } else if ( at == TOUCH_21_INTERIOR_2 ) {
	return TOUCH_12_INTERIOR_1;
      } else if ( at == TOUCH_22_INTERIOR_1 ) {
	return TOUCH_22_INTERIOR_2;
      } else if ( at == TOUCH_22_INTERIOR_2 ) {
	return TOUCH_22_INTERIOR_1;
      } else if ( at == OVERLAPPING_12 ) {
	return OVERLAPPING_21;
      } else if ( at == OVERLAPPING_21 ) {
	return OVERLAPPING_12;
      } else if ( at == TOUCH_INTERIOR_12 ) {
	return TOUCH_INTERIOR_11;
      } else if ( at == TOUCH_INTERIOR_22 ) {
	return TOUCH_INTERIOR_21;
      } else if ( at == TOUCH_INTERIOR_11 ) {
	return TOUCH_INTERIOR_12;
      } else if ( at == TOUCH_INTERIOR_21 ) {
	return TOUCH_INTERIOR_22;
      }
      return at;
    }
  };

#ifdef CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARRANGEMENT_TYPE_OUTPUT_OPERATOR
  static
  std::ostream& operator<<(std::ostream& os,
			   const Arrangement_enum::Arrangement_type& at)
  {
    typedef Arrangement_enum AT;

    if ( at == AT::DISJOINT ) {
      os << "DISJOINT";
    } else if ( at == AT::TOUCH_1 ) {
      os << "TOUCH_1";
    } else if ( at == AT::TOUCH_2 ) {
      os << "TOUCH_2";
    } else if ( at == AT::TOUCH_11 ) {
      os << "TOUCH_11";
    } else if ( at == AT::TOUCH_12 ) {
      os << "TOUCH_12";
    } else if ( at == AT::TOUCH_21 ) {
      os << "TOUCH_21";
    } else if ( at == AT::TOUCH_22 ) {
      os << "TOUCH_22";
    } else if ( at == AT::CROSSING ) {
      os << "CROSSING";
    } else if ( at == AT::IDENTICAL) {
      os << "IDENTICAL";
    } else if ( at == AT::INTERIOR_1 ) {
      os << "INTERIOR_1";
    } else if ( at == AT::INTERIOR_2 ) {
      os << "INTERIOR_2";
    } else if ( at == AT::INTERIOR ) {
      os << "INTERIOR";
    } else if ( at == AT::TOUCH_11_INTERIOR_1 ) {
      os << "TOUCH_11_INTERIOR_1";
    } else if ( at == AT::TOUCH_11_INTERIOR_2 ) {
      os << "TOUCH_11_INTERIOR_2";
    } else if ( at == AT::TOUCH_12_INTERIOR_1 ) {
      os << "TOUCH_12_INTERIOR_1";
    } else if ( at == AT::TOUCH_12_INTERIOR_2 ) {
      os << "TOUCH_12_INTERIOR_2";
    } else if ( at == AT::TOUCH_21_INTERIOR_1 ) {
      os << "TOUCH_21_INTERIOR_1";
    } else if ( at == AT::TOUCH_21_INTERIOR_2 ) {
      os << "TOUCH_21_INTERIOR_2";
    } else if ( at == AT::TOUCH_22_INTERIOR_1 ) {
      os << "TOUCH_22_INTERIOR_1";
    } else if ( at == AT::TOUCH_22_INTERIOR_2 ) {
      os << "TOUCH_22_INTERIOR_2";
    } else if ( at == AT::OVERLAPPING_11 ) {
      os << "OVERLAPPING_11";
    } else if ( at == AT::OVERLAPPING_12 ) {
      os << "OVERLAPPING_12";
    } else if ( at == AT::OVERLAPPING_21 ) {
      os << "OVERLAPPING_21";
    } else if ( at == AT::OVERLAPPING_22 ) {
      os << "OVERLAPPING_22";
    } else if ( at == AT::TOUCH_INTERIOR_11 ) {
      os << "TOUCH_INTERIOR_11";
    } else if ( at == AT::TOUCH_INTERIOR_12 ) {
      os << "TOUCH_INTERIOR_12";
    } else if ( at == AT::TOUCH_INTERIOR_21 ) {
      os << "TOUCH_INTERIOR_21";
    } else if ( at == AT::TOUCH_INTERIOR_22 ) {
      os << "TOUCH_INTERIOR_22";
    } else {
      CGAL_error();
    }

    return os;
  }
#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARRANGEMENT_TYPE_OUTPUT_OPERATOR

} // namespace Internal


} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_ARRANGEMENT_ENUM_H
