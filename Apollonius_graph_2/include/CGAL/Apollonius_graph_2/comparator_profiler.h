// Copyright (c) 2003,2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_COMPARATOR_PROFILER_H
#define CGAL_COMPARATOR_PROFILER_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Apollonius_graph_2/basic.h>
#include <CGAL/atomic.h>

namespace CGAL {

namespace ApolloniusGraph_2 {

class comparator_profiler
{
public:

#ifdef CGAL_NO_ATOMIC
  typedef bool bool_;
  typedef unsigned long long_;
#else
  typedef CGAL::cpp11::atomic<bool> bool_;
  typedef CGAL::cpp11::atomic<unsigned long> long_;
#endif

  static bool_ count_cases;
  static long_ case_1_counter;
  static long_ case_2_counter;
  static long_ case_3a_Jpos_counter;
  static long_ case_3a_Jneg_counter;
  static long_ case_3b_Jpos_counter;
  static long_ case_3b_Jneg_counter;
  static long_ case_4_counter;
  static long_ case_5_counter;
  static long_ case_degenerate_counter;
public:
  static long_ counter_rr;
  static long_ counter_rr_p3inf;
  static long_ counter_rr_p4;
  static long_ counter_rr_e;
  static long_ counter_rr_r0;

  static void reset()
  {
    count_cases = false;
    case_1_counter = 0;
    case_2_counter = 0;
    case_3a_Jpos_counter = 0;
    case_3a_Jneg_counter = 0;
    case_3b_Jpos_counter = 0;
    case_3b_Jneg_counter = 0;
    case_4_counter = 0;
    case_5_counter = 0;
    case_degenerate_counter = 0;

    counter_rr = 0;
    counter_rr_p3inf = 0;
    counter_rr_p4 = 0;
    counter_rr_e = 0;
    counter_rr_r0 = 0;
  }

  template< class FT >
  static void count_case(const FT& a1, const FT& b1, const FT& c1,
			 const FT& a2, const FT& b2, const FT& c2)
  {
    // works correctly only with leda_real
    FT D1 = CGAL::square(b1) - a1 * c1;
    
    FT l1 = (b1 - CGAL::sqrt(D1)) / a1;
    FT r1 = (b1 + CGAL::sqrt(D1)) / a1;
    if ( a1 < 0 ) { std::swap(r1, l1); }
    
    FT D2 = CGAL::square(b2) - a2 * c2;

    if ( D1 == 0 || D2 == 0 ) {
      case_degenerate_counter++;
      return;
    }

    FT l2 = (b2 - CGAL::sqrt(D2)) / a2;
    FT r2 = (b2 + CGAL::sqrt(D2)) / a2;
    if ( a2 < 0 ) { std::swap(r2, l2); }

    if ( l1 < l2 ) {
      if ( r1 > r2 ) {
	FT J = a1 * b2 - a2 * b1;
	if ( J > 0 ) {
	  case_3b_Jpos_counter++;
	} else if ( J < 0 ) {
	  case_3b_Jneg_counter++;
	} else {
	  case_degenerate_counter++;
	}
      } else if ( r1 < r2 ) {
	if ( r1 < l2 ) {
	  case_5_counter++;
	} else if ( r1 > l2 ) {
	  case_4_counter++;
	} else {
	  case_degenerate_counter++;
	}
      } else {
	case_degenerate_counter++;
      }
    } else if ( l1 > l2 ) {
      if ( r1 < r2 ) {
	FT J = a1 * b2 - a2 * b1;
	if ( J > 0 ) {
	  case_3a_Jpos_counter++;
	} else if ( J < 0 ) {
	  case_3a_Jneg_counter++;
	} else {
	  case_degenerate_counter++;
	}
      } else if ( r1 > r2 ) {
	if ( l1 < r2 ) {
	  case_2_counter++;
	} else if ( l1 > r2 ) {
	  case_1_counter++;
	} else {
	  case_degenerate_counter++;
	}
      } else {
	case_degenerate_counter++;
      }
    } else {
      case_degenerate_counter++;
    }
  }
};

#ifdef CGAL_NO_ATOMIC

bool comparator_profiler::count_cases = false;
unsigned long comparator_profiler::case_1_counter = 0;
unsigned long comparator_profiler::case_2_counter = 0;
unsigned long comparator_profiler::case_3a_Jpos_counter = 0;
unsigned long comparator_profiler::case_3a_Jneg_counter = 0;
unsigned long comparator_profiler::case_3b_Jpos_counter = 0;
unsigned long comparator_profiler::case_3b_Jneg_counter = 0;
unsigned long comparator_profiler::case_4_counter = 0;
unsigned long comparator_profiler::case_5_counter = 0;
unsigned long comparator_profiler::case_degenerate_counter = 0;

unsigned long comparator_profiler::counter_rr = 0;
unsigned long comparator_profiler::counter_rr_p3inf = 0;
unsigned long comparator_profiler::counter_rr_p4 = 0;
unsigned long comparator_profiler::counter_rr_e = 0;
unsigned long comparator_profiler::counter_rr_r0 = 0;

#endif


} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_COMPARATOR_PROFILER_H
