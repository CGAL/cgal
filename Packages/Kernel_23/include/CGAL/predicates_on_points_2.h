// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : predicates_on_points_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_PREDICATES_ON_POINTS_2_H
#define CGAL_PREDICATES_ON_POINTS_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/predicates_on_pointsH2.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/predicates_on_points_2.h>
#endif

#include <CGAL/Point_2.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
bool x_equal(const Point_2<R>& p,
             const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
  return x_equal(static_cast<const RPoint_2&>(p),
	         static_cast<const RPoint_2&>(q));
}

template < class R >
inline
bool y_equal(const Point_2<R>& p,
             const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
  return y_equal(static_cast<const RPoint_2&>(p),
	         static_cast<const RPoint_2&>(q));
}


template < class R >
inline
Comparison_result compare_x(const Point_2<R>& p,
                            const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
  return compare_x(static_cast<const RPoint_2&>(p),
	           static_cast<const RPoint_2&>(q));
}

template < class R >
inline
Comparison_result compare_y(const Point_2<R>& p,
                            const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
  return compare_y(static_cast<const RPoint_2&>(p),
	           static_cast<const RPoint_2&>(q));
}

template < class R >
inline
Comparison_result
compare_deltax_deltay(const Point_2<R>& p,
                      const Point_2<R>& q,
                      const Point_2<R>& r,
                      const Point_2<R>& s)
{
  typedef typename R::Point_2_base  RPoint_2;
  return compare_deltax_deltay(static_cast<const RPoint_2&>(p),
                               static_cast<const RPoint_2&>(q),
                               static_cast<const RPoint_2&>(r),
                               static_cast<const RPoint_2&>(s));
}

template < class R >
inline
Comparison_result
compare_lexicographically_xy(const Point_2<R>& p,
                             const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
  return  compare_lexicographically_xy(static_cast<const RPoint_2&>(p),
                                       static_cast<const RPoint_2&>(q));
}

template < class R >
inline
bool
lexicographically_xy_smaller_or_equal(const Point_2<R>& p,
                                      const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
  return lexicographically_xy_smaller_or_equal(
	         static_cast<const RPoint_2&>(p),
		 static_cast<const RPoint_2&>(q));
}

template < class R >
inline
bool
lexicographically_xy_smaller(const Point_2<R>& p,
                                  const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
 return lexicographically_xy_smaller(static_cast<const RPoint_2&>(p),
                                     static_cast<const RPoint_2&>(q));
}

template < class R >
inline
bool
lexicographically_xy_larger_or_equal(const Point_2<R>& p,
                                     const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
 return !lexicographically_xy_smaller(static_cast<const RPoint_2&>(p),
                                      static_cast<const RPoint_2&>(q));
}

template < class R >
inline
bool
lexicographically_xy_larger(const Point_2<R>& p,
                            const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
  return !lexicographically_xy_smaller_or_equal(
	               static_cast<const RPoint_2&>(p),
		       static_cast<const RPoint_2&>(q));
}

template < class R >
inline
Comparison_result
compare_lexicographically_yx(const Point_2<R>& p,
                                  const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
  return  compare_lexicographically_yx(static_cast<const RPoint_2&>(p),
                                       static_cast<const RPoint_2&>(q));
}

template < class R >
inline
bool
lexicographically_yx_smaller_or_equal(const Point_2<R>& p,
                                      const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
  return lexicographically_yx_smaller_or_equal(
	            static_cast<const RPoint_2&>(p),
		    static_cast<const RPoint_2&>(q));
}

template < class R >
inline
bool
lexicographically_yx_smaller(const Point_2<R>& p,
                             const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
 return lexicographically_yx_smaller(static_cast<const RPoint_2&>(p),
                                     static_cast<const RPoint_2&>(q));
}

template < class R >
inline
bool
lexicographically_yx_larger_or_equal(const Point_2<R>& p,
                                     const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
 return !lexicographically_yx_smaller(static_cast<const RPoint_2&>(p),
                                      static_cast<const RPoint_2&>(q));
}

template < class R >
inline
bool
lexicographically_yx_larger(const Point_2<R>& p,
                            const Point_2<R>& q)
{
  typedef typename R::Point_2_base  RPoint_2;
  return !lexicographically_yx_smaller_or_equal(
	          static_cast<const RPoint_2&>(p),
		  static_cast<const RPoint_2&>(q));
}

template < class R >
inline
bool
are_ordered_along_line(const Point_2<R>& p,
                       const Point_2<R>& q,
                       const Point_2<R>& r)
{
  typedef typename R::Point_2_base  RPoint_2;
  return are_ordered_along_line(static_cast<const RPoint_2&>(p),
                                static_cast<const RPoint_2&>(q),
                                static_cast<const RPoint_2&>(r));
}

template < class R >
inline
bool
collinear_are_ordered_along_line(const Point_2<R>& p,
                                 const Point_2<R>& q,
                                 const Point_2<R>& r)
{
  typedef typename R::Point_2_base  RPoint_2;
  return collinear_are_ordered_along_line(static_cast<const RPoint_2&>(p),
                                          static_cast<const RPoint_2&>(q),
                                          static_cast<const RPoint_2&>(r));
}

template < class R >
inline
bool
are_strictly_ordered_along_line(const Point_2<R>& p,
                            const Point_2<R>& q,
                            const Point_2<R>& r)
{
  typedef typename R::Point_2_base  RPoint_2;
  return are_strictly_ordered_along_line(static_cast<const RPoint_2&>(p),
                                         static_cast<const RPoint_2&>(q),
                                         static_cast<const RPoint_2&>(r));
}
template < class R >
inline
bool
collinear_are_strictly_ordered_along_line(const Point_2<R>& p,
                                          const Point_2<R>& q,
                                          const Point_2<R>& r)
{
  typedef typename R::Point_2_base  RPoint_2;
  return
  collinear_are_strictly_ordered_along_line(static_cast<const RPoint_2&>(p),
                                            static_cast<const RPoint_2&>(q),
                                            static_cast<const RPoint_2&>(r));
}
template < class R >
inline
bool
collinear(const Point_2<R>& p,
          const Point_2<R>& q,
          const Point_2<R>& r)
{
  typedef typename R::Point_2_base  RPoint_2;
  return (collinear(static_cast<const RPoint_2&>(p),
                    static_cast<const RPoint_2&>(q),
                    static_cast<const RPoint_2&>(r)));
}

template < class R >
inline
bool
leftturn(const Point_2<R>& p,
         const Point_2<R>& q,
         const Point_2<R>& r)
{
  typedef typename R::Point_2_base  RPoint_2;
  return leftturn(static_cast<const RPoint_2&>(p),
                  static_cast<const RPoint_2&>(q),
                  static_cast<const RPoint_2&>(r));
}

template < class R >
inline
bool
rightturn(const Point_2<R>& p,
          const Point_2<R>& q,
          const Point_2<R>& r)
{
  typedef typename R::Point_2_base  RPoint_2;
  return rightturn(static_cast<const RPoint_2&>(p),
                   static_cast<const RPoint_2&>(q),
                   static_cast<const RPoint_2&>(r));
}

template < class R >
inline
bool
rightturn(const Origin& o,
          const Point_2<R>& q,
          const Point_2<R>& r)
{
  typedef typename R::Point_2_base  RPoint_2;
  return rightturn(o, static_cast<const RPoint_2&>(q),
	              static_cast<const RPoint_2&>(r));
}

template < class R >
inline
Orientation
orientation(const Point_2<R>& p,
            const Point_2<R>&q,
            const Point_2<R>& r)
{
  typedef typename R::Point_2_base  RPoint_2;
  return orientation(static_cast<const RPoint_2&>(p),
                     static_cast<const RPoint_2&>(q),
                     static_cast<const RPoint_2&>(r));
}

template <class R >
inline
Oriented_side
side_of_oriented_circle(const Point_2<R>& p,
                        const Point_2<R>& q,
                        const Point_2<R>& r,
                        const Point_2<R>& test)
{
  typedef typename R::Point_2_base  RPoint_2;
  return side_of_oriented_circle(static_cast<const RPoint_2&>(p),
                                 static_cast<const RPoint_2&>(q),
                                 static_cast<const RPoint_2&>(r),
                                 static_cast<const RPoint_2&>(test));
}

template <class R >
inline
Bounded_side
side_of_bounded_circle(const Point_2<R>& p,
                       const Point_2<R>& q,
                       const Point_2<R>& r,
                       const Point_2<R>& test)
{
  typedef typename R::Point_2_base  RPoint_2;
  return side_of_bounded_circle(static_cast<const RPoint_2&>(p),
                                static_cast<const RPoint_2&>(q),
                                static_cast<const RPoint_2&>(r),
                                static_cast<const RPoint_2&>(test));
}

CGAL_END_NAMESPACE

#endif  // CGAL_PREDICATES_ON_POINTS_2_H
