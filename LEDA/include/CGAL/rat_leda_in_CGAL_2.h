// Copyright (c) 1999  
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
// 
//
// Author(s)     : Stefan Schirra


#ifndef CGAL_RAT_LEDA_IN_CGAL_2_H
#define CGAL_RAT_LEDA_IN_CGAL_2_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/LEDA_basic.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/rat_point.h>
#include <LEDA/rat_segment.h>
#include <LEDA/rat_line.h>
#else
#include <LEDA/geo/rat_point.h>
#include <LEDA/geo/rat_segment.h>
#include <LEDA/geo/rat_line.h>
#endif

namespace CGAL {

inline
bool
left_turn( const leda_rat_point & p, 
          const leda_rat_point & q, 
          const leda_rat_point & r)
{ return CGAL_LEDA_SCOPE::left_turn(p,q,r); }

inline
bool
right_turn( const leda_rat_point & p, 
                const leda_rat_point & q, 
                const leda_rat_point & r)
{ return CGAL_LEDA_SCOPE::right_turn(p,q,r); }

/*
inline
bool
collinear( const leda_rat_point & p, 
           const leda_rat_point & q, 
           const leda_rat_point & r)
{ return ::collinear(p,q,r); }
*/


inline
Orientation
orientation( const leda_rat_point & p, 
             const leda_rat_point & q, 
             const leda_rat_point & r)
{ return (Orientation)CGAL_LEDA_SCOPE::orientation(p,q,r); }


inline
bool
lexicographically_xy_smaller( const leda_rat_point & p, 
                              const leda_rat_point & q)
{ return ( leda_rat_point::cmp_xy(p,q)  <  0 ); }

inline
bool
lexicographically_yx_smaller( const leda_rat_point & p, 
                              const leda_rat_point & q)
{ return ( leda_rat_point::cmp_yx(p,q)  <  0 ); }

inline
bool
lexicographically_xy_larger( const leda_rat_point & p, 
                             const leda_rat_point & q)
{ return ( leda_rat_point::cmp_xy(p,q)  >  0 ); }

inline
bool
lexicographically_yx_larger( const leda_rat_point & p, 
                             const leda_rat_point & q)
{ return ( leda_rat_point::cmp_yx(p,q)  >  0 ); }

inline
bool
lexicographically_xy_smaller_or_equal( const leda_rat_point & p, 
                                       const leda_rat_point & q)
{ return ( leda_rat_point::cmp_xy(p,q)  <=  0 ); }

inline
bool
lexicographically_yx_smaller_or_equal( const leda_rat_point & p, 
                                       const leda_rat_point & q)
{ return ( leda_rat_point::cmp_yx(p,q)  <=  0 ); }

inline 
Comparison_result
compare_lexicographically_xy(const leda_rat_point & p, 
                             const leda_rat_point & q)
{ return (Comparison_result)leda_rat_point::cmp_xy(p,q); }

inline 
Comparison_result
compare_lexicographically_yx(const leda_rat_point & p, 
                             const leda_rat_point & q)
{ return (Comparison_result)leda_rat_point::cmp_yx(p,q); }

inline
bool
collinear_are_ordered_along_line( const leda_rat_point & p, 
                                  const leda_rat_point & q, 
                                  const leda_rat_point & r)
{ return leda_rat_segment(p,r).contains(q); }  

inline
bool
collinear_are_strictly_ordered_along_line( const leda_rat_point & p, 
                                           const leda_rat_point & q, 
                                           const leda_rat_point & r)
{ return (leda_rat_segment(p,r).contains(q) && ( q != p ) && ( q != r )); }  

inline
bool
are_ordered_along_line( const leda_rat_point & p, 
                        const leda_rat_point & q, 
                        const leda_rat_point & r)
{ return ( CGAL_LEDA_SCOPE::collinear(p,q,r) &&  collinear_are_ordered_along_line(p,q,r)); }  

inline
bool
are_strictly_ordered_along_line( const leda_rat_point & p, 
                                 const leda_rat_point & q, 
                                 const leda_rat_point & r)
{ 
  return (    CGAL_LEDA_SCOPE::collinear(p,q,r) 
           && collinear_are_strictly_ordered_along_line(p,q,r)
         ); 
}  

inline
Comparison_result 
cmp_signed_dist_to_line(const leda_rat_point& p, const leda_rat_point& q,
                        const leda_rat_point& r, const leda_rat_point& s)
{
#if ( __LEDA__ >= 360 )
  return (Comparison_result)CGAL_LEDA_SCOPE::cmp_signed_dist(p,q,r,s);
#else
  leda_rat_line  l(p,q);
  int  r_or = CGAL_LEDA_SCOPE::orientation( l, r );
  int  s_or = CGAL_LEDA_SCOPE::orientation( l, s );
  if ( r_or != s_or )
  {
      return (Comparison_result)( r_or < s_or );
  }
  else
  {
     return 
      (Comparison_result)(r_or *( CGAL::sign(l.sqr_dist(r) - l.sqr_dist(s) )));
  }
#endif  // __LEDA__ >= 360
}

inline
bool
has_smaller_signed_dist_to_line(const leda_rat_point& p, 
                                const leda_rat_point& q,
                                const leda_rat_point& r, 
                                const leda_rat_point& s)
{ return ( cmp_signed_dist_to_line(p,q,r,s) == SMALLER ); }

inline
bool
has_larger_signed_dist_to_line(const leda_rat_point& p, 
                               const leda_rat_point& q,
                               const leda_rat_point& r, 
                               const leda_rat_point& s)
{ return ( cmp_signed_dist_to_line(p,q,r,s) == LARGER ); }

} //namespace CGAL

#endif // CGAL_RAT_LEDA_IN_CGAL_2_H

