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
// file          : include/CGAL/convex_hull_rat_leda_traits_2.h
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

#ifndef CGAL_CONVEX_HULL_RAT_LEDA_TRAITS_2_H
#define CGAL_CONVEX_HULL_RAT_LEDA_TRAITS_2_H

#ifndef CGAL_NO_DEPRECATED_CODE

#include <CGAL/basic.h>
#include <CGAL/rat_leda_in_CGAL_2.h>
#include <CGAL/Kernel/function_objects.h>

CGAL_BEGIN_NAMESPACE

struct LEDA_equal_rat_points {
  typedef Arity_tag< 2 > Arity;
  typedef bool           result_type;  

  bool operator()(const leda_rat_point& p1, const leda_rat_point& p2) const
  { return (p1 == p2); }
};

class LEDA_rat_kernel_2
{
public:
   typedef leda_rat_point     Point_2;
   typedef leda_rat_segment   Segment_2;
};

class Convex_hull_rat_leda_traits_2
{
public:
  typedef   leda_rat_point                                  Point_2;    
  typedef   CGALi::Less_xy_2<LEDA_rat_kernel_2>             Less_xy_2;
  typedef   CGALi::Less_yx_2<LEDA_rat_kernel_2>             Less_yx_2;
  typedef   CGALi::Less_signed_distance_to_line_2<LEDA_rat_kernel_2>
                                               Less_signed_distance_to_line_2;
  typedef   CGALi::Less_rotate_ccw_2<LEDA_rat_kernel_2>     Less_rotate_ccw_2;
  typedef   CGALi::Left_turn_2<LEDA_rat_kernel_2>           Left_turn_2;
  typedef   CGAL::LEDA_equal_rat_points                     Equal_2;  
  typedef   leda_rat_segment                                Segment_2; 
  
  Less_xy_2
  less_xy_2_object() const 
  { return Less_xy_2(); } 

  Less_yx_2
  less_yx_2_object() const 
  { return Less_yx_2(); } 

  Less_signed_distance_to_line_2
  less_signed_distance_to_line_2_object( ) const
  { return Less_signed_distance_to_line_2( ); } 

  Less_rotate_ccw_2
  less_rotate_ccw_2_object( ) const
  { return Less_rotate_ccw_2( ); }

  Left_turn_2
  left_turn_2_object() const
  { return Left_turn_2(); }

  Equal_2
  equal_2_object() const
  { return Equal_2(); }
};

// for backward compatability
typedef Convex_hull_rat_leda_traits_2   convex_hull_rat_leda_traits_2;

CGAL_END_NAMESPACE

#endif // CGAL_NO_DEPRECATED_CODE

#endif // CGAL_cONVEX_HULL_RAT_LEDA_TRAITS_2_H

