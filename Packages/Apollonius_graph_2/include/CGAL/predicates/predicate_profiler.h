// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/predicates/predicate_profiler.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



#ifndef CGAL_AG2_PREDICATE_PROFILER
#define CGAL_AG2_PREDICATE_PROFILER

#define AG2_PROFILE_PREDICATES

CGAL_BEGIN_NAMESPACE

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

CGAL_END_NAMESPACE

#endif // CGAL_AG2_PREDICATE_PROFILER
