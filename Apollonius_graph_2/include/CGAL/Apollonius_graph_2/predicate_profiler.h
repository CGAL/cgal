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



#ifndef CGAL_APOLLONIUS_GRAPH_2_PREDICATE_PROFILER
#define CGAL_APOLLONIUS_GRAPH_2_PREDICATE_PROFILER 1

#include <CGAL/Apollonius_graph_2/basic.h>

#define AG2_PROFILE_PREDICATES

namespace CGAL {

namespace ApolloniusGraph_2 {

class ag2_predicate_profiler
{
public:
  // high level predicates
  static unsigned long side_of_bisector_counter;
  static unsigned long is_trivial_counter;
  static unsigned long infinite_edge_conflict_type_counter;
  static unsigned long finite_edge_conflict_type_counter;

  // subpredicates
  static unsigned long inside_circular_arc_counter;
  static unsigned long distance_from_bitangent_counter;
  static unsigned long shadow_region_type_counter;
  static unsigned long incircle_counter;
  static unsigned long order_on_bisector_counter;

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

unsigned long ag2_predicate_profiler::side_of_bisector_counter = 0;
unsigned long ag2_predicate_profiler::is_trivial_counter = 0;
unsigned long ag2_predicate_profiler::infinite_edge_conflict_type_counter = 0;
unsigned long ag2_predicate_profiler::finite_edge_conflict_type_counter = 0;

unsigned long ag2_predicate_profiler::inside_circular_arc_counter = 0;
unsigned long ag2_predicate_profiler::distance_from_bitangent_counter = 0;
unsigned long ag2_predicate_profiler::shadow_region_type_counter = 0;
unsigned long ag2_predicate_profiler::incircle_counter = 0;
unsigned long ag2_predicate_profiler::order_on_bisector_counter = 0;

} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_PREDICATE_PROFILER
