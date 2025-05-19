// Copyright (c) 2003,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_APOLLONIUS_GRAPH_2_PREDICATE_PROFILER
#define CGAL_APOLLONIUS_GRAPH_2_PREDICATE_PROFILER 1

#include <CGAL/license/Apollonius_graph_2.h>


#include <CGAL/Apollonius_graph_2/basic.h>
#include <atomic>

#define AG2_PROFILE_PREDICATES

namespace CGAL {

namespace ApolloniusGraph_2 {

class ag2_predicate_profiler
{
public:
#ifdef CGAL_NO_ATOMIC
  typedef unsigned long long_;
#else
  typedef std::atomic<unsigned long> long_;
#endif

  // high level predicates
  static long_ side_of_bisector_counter;
  static long_ is_trivial_counter;
  static long_ infinite_edge_conflict_type_counter;
  static long_ finite_edge_conflict_type_counter;

  // subpredicates
  static long_ inside_circular_arc_counter;
  static long_ distance_from_bitangent_counter;
  static long_ shadow_region_type_counter;
  static long_ incircle_counter;
  static long_ order_on_bisector_counter;

  static void reset() {
    side_of_bisector_counter = 0;
    is_trivial_counter = 0;
    infinite_edge_conflict_type_counter = 0;
    finite_edge_conflict_type_counter = 0;

    inside_circular_arc_counter = 0;
    distance_from_bitangent_counter = 0;
    shadow_region_type_counter = 0;
    incircle_counter = 0;
    order_on_bisector_counter = 0;
  }
};


#ifdef CGAL_NO_ATOMIC

unsigned long ag2_predicate_profiler::side_of_bisector_counter = 0;
unsigned long ag2_predicate_profiler::is_trivial_counter = 0;
unsigned long ag2_predicate_profiler::infinite_edge_conflict_type_counter = 0;
unsigned long ag2_predicate_profiler::finite_edge_conflict_type_counter = 0;

unsigned long ag2_predicate_profiler::inside_circular_arc_counter = 0;
unsigned long ag2_predicate_profiler::distance_from_bitangent_counter = 0;
unsigned long ag2_predicate_profiler::shadow_region_type_counter = 0;
unsigned long ag2_predicate_profiler::incircle_counter = 0;
unsigned long ag2_predicate_profiler::order_on_bisector_counter = 0;

#endif

} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_PREDICATE_PROFILER
