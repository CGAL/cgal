// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : include/CGAL/convex_hull_constructive_traits_2.h
// package       : Convex_hull_2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

// This file's name must begin with a lower-case letter for backward 
// compatability.  Unfortunately, you can't have a file that differs only 
// in capitalization on the Windows platforms.

#ifndef CGAL_CONVEX_HULL_CONSTRUCTIVE_TRAITS_2_H 
#define CGAL_CONVEX_HULL_CONSTRUCTIVE_TRAITS_2_H


#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/distance_predicates_2.h>
#include <CGAL/predicate_objects_on_points_2.h>

CGAL_BEGIN_NAMESPACE
template <class K_>
class Convex_hull_constructive_traits_2 : public K_
{
public:
  typedef K_                                  K;
  typedef typename K::Point_2                 Point_2;    
  typedef typename K::Less_xy_2               Less_xy_2;
  typedef typename K::Less_yx_2               Less_yx_2;
  typedef CGAL::r_Less_dist_to_line<K>        Less_signed_distance_to_line_2;
  typedef typename K::Less_rotate_ccw_2       Less_rotate_ccw_2;
  typedef typename K::Leftturn_2              Leftturn_2;
  typedef typename K::Segment_2               Segment_2;    
  
  Less_xy_2
  less_xy_2_object() const 
  { return Less_xy_2(); } 

  Less_yx_2
  less_yx_2_object() const 
  { return Less_yx_2(); } 

  Less_signed_distance_to_line_2
  less_signed_distance_to_line_2_object() const
  { return Less_signed_distance_to_line_2(); } 

  Less_rotate_ccw_2
  less_rotate_ccw_2_object() const
  { return Less_rotate_ccw_2(); }

  Leftturn_2
  leftturn_2_object() const
  { return Leftturn_2(); }

};


// for backward compatability

template <class K>
class convex_hull_constructive_traits_2 : 
                         public Convex_hull_constructive_traits_2<K>
{
};

CGAL_END_NAMESPACE

#endif // CGAL_CONVEX_HULL_CONSTRUCTIVE_TRAITS_2_H

