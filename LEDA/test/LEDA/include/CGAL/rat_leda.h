// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_RAT_LEDA_H
#define CGAL_RAT_LEDA_H

#include <CGAL/LEDA_basic.h>
#include <CGAL/user_classes.h>
#include <CGAL/basic_classes.h>
#if CGAL_LEDA_VERSION < 500
#include <LEDA/rat_point.h>
#else
#include <LEDA/geo/rat_point.h>
#endif
#include <CGAL/Homogeneous.h>

namespace CGAL {

class use_rat_leda_kernel;
template <>
class Point_2<use_rat_leda_kernel>;

class use_rat_leda_kernel : public Homogeneous_base< use_rat_leda_kernel, leda_integer, leda_rational>
{
  public:
    typedef  use_rat_leda_kernel  R;

    typedef  leda_integer      RT;
    typedef  leda_rational     FT;

    typedef  leda_rat_point    Point_2_base;
    typedef  leda_rat_vector   Vector_2_base;

    typedef  CGAL::Point_2<R>  Point_2;
    typedef  CGAL::Vector_2<R> Vector_2;
};


template <>
class Point_2<use_rat_leda_kernel> : public leda_rat_point
{
  public:

    typedef use_rat_leda_kernel      R;
    typedef use_rat_leda_kernel::RT  RT;
    typedef use_rat_leda_kernel::FT  FT;

    Point_2<use_rat_leda_kernel>()
    {}

    Point_2<use_rat_leda_kernel>(const Origin &o)
    {}

    Point_2<use_rat_leda_kernel>(const leda_rat_point& p)
      : leda_rat_point(p)
    {}

    Point_2<use_rat_leda_kernel>(const RT &hx, const RT &hy)
      : leda_rat_point(hx, hy)
    {}

    Point_2<use_rat_leda_kernel>(const RT &hx, const RT &hy, const RT &hw)
      : leda_rat_point(hx, hy, hw)
    {}
};

inline
Comparison_result
compare_distance_to_point( const leda_rat_point& p,
                   const leda_rat_point& q,
                   const leda_rat_point& r)
{
  if ( q == r) return EQUAL;
  return ( p.sqr_dist(q) > p.sqr_dist(r) ) ? LARGER : SMALLER;
}

inline
bool
has_larger_distance_to_point( const leda_rat_point& p,
                          const leda_rat_point& q,
                          const leda_rat_point& r)
{ return (compare_distance_to_point(p,q,r) == LARGER); }

inline
bool
has_smaller_distance_to_point( const leda_rat_point& p,
                           const leda_rat_point& q,
                           const leda_rat_point& r)
{ return (compare_distance_to_point(p,q,r) == SMALLER); }

inline
bool
has_smaller_signed_distance_to_line( const leda_rat_point& p,
                                     const leda_rat_point& q,
                                     const leda_rat_point& r,
                                     const leda_rat_point& s)
{ return (cmp_signed_dist(p,q,r,s) == SMALLER); }

inline
bool
has_larger_signed_distance_to_line( const leda_rat_point& p,
                                     const leda_rat_point& q,
                                     const leda_rat_point& r,
                                     const leda_rat_point& s)
{ return (cmp_signed_dist(p,q,r,s) == LARGER); }

inline
Comparison_result
compare_signed_distance_to_line( const leda_rat_point& p,
                                     const leda_rat_point& q,
                                     const leda_rat_point& r,
                                     const leda_rat_point& s)
{ return Comparison_result(cmp_signed_dist(p,q,r,s)); }

inline
bool 
left_turn(const leda_rat_point& p, const leda_rat_point& q, 
          const leda_rat_point& r)
{
   return CGAL_LEDA_SCOPE::left_turn(p,q,r);
}

inline
bool 
right_turn(const leda_rat_point& p, const leda_rat_point& q, 
          const leda_rat_point& r)
{
   return CGAL_LEDA_SCOPE::right_turn(p,q,r);
}

inline
bool
lexicographically_yx_larger_or_equal(const leda_rat_point& p,
                                     const leda_rat_point& q)
{
   if (p.ycoord() > q.ycoord()) return true;
   if (p.ycoord() < q.ycoord()) return false;
   return (p.xcoord() >= q.xcoord()); 
}

inline
bool
lexicographically_xy_larger_or_equal(const leda_rat_point& p,
                                     const leda_rat_point& q)
{
   if (p.xcoord() > q.xcoord()) return true;
   if (p.xcoord() < q.xcoord()) return false;
   return (p.ycoord() >= q.ycoord()); 
}

inline
Comparison_result
compare_yx(const leda_rat_point& p, const leda_rat_point& q)
{
   if (p.ycoord() > q.ycoord()) return CGAL::LARGER;
   if (p.ycoord() < q.ycoord()) return CGAL::SMALLER;
   if (p.xcoord() > q.xcoord()) return CGAL::LARGER; 
   if (p.xcoord() < q.xcoord()) return CGAL::SMALLER; 
   return CGAL::EQUAL;
}

inline
Comparison_result
compare_xy(const leda_rat_point& p, const leda_rat_point& q)
{
   if (p.xcoord() > q.xcoord()) return CGAL::LARGER; 
   if (p.xcoord() < q.xcoord()) return CGAL::SMALLER; 
   if (p.ycoord() > q.ycoord()) return CGAL::LARGER;
   if (p.ycoord() < q.ycoord()) return CGAL::SMALLER;
   return CGAL::EQUAL;
}
} //namespace CGAL

#endif // CGAL_RAT_LEDA_H
