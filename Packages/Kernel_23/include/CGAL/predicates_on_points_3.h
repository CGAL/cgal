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
// file          : predicates_on_points_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_PREDICATES_ON_POINTS_3_H
#define CGAL_PREDICATES_ON_POINTS_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/predicates_on_pointsH3.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/predicates_on_points_3.h>
#endif

#include <CGAL/Point_3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Comparison_result
compare_lexicographically_xyz( const Point_3<R> &p,
                               const Point_3<R> &q)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return compare_lexicographically_xyz(static_cast<const RPoint_3&>(p),
                                       static_cast<const RPoint_3&>(q));
}

template < class R >
inline
bool
lexicographically_xyz_smaller_or_equal(const Point_3<R> &p,
                                       const Point_3<R> &q)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return ( !( compare_lexicographically_xyz(static_cast<const RPoint_3&>(p),
                                            static_cast<const RPoint_3&>(q))
              == LARGER ) );
}

template < class R >
inline
bool
lexicographically_xyz_smaller(const Point_3<R> &p,
                              const Point_3<R> &q)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return ( compare_lexicographically_xyz(static_cast<const RPoint_3&>(p),
                                         static_cast<const RPoint_3&>(q))
           == SMALLER );
}

template < class R >
inline
Comparison_result
compare_x(const Point_3<R> &p, const Point_3<R> &q)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return compare_x(static_cast<const RPoint_3&>(p),
	           static_cast<const RPoint_3&>(q));
}


template < class R >
inline
Comparison_result
compare_y(const Point_3<R> &p, const Point_3<R> &q)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return compare_y(static_cast<const RPoint_3&>(p),
	           static_cast<const RPoint_3&>(q));
}


template < class R >
inline
Comparison_result
compare_z(const Point_3<R> &p, const Point_3<R> &q)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return compare_z(static_cast<const RPoint_3&>(p),
	           static_cast<const RPoint_3&>(q));
}

template < class R >
inline
bool
x_equal(const Point_3<R> &p,
        const Point_3<R> &q)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return x_equal(static_cast<const RPoint_3&>(p),
	         static_cast<const RPoint_3&>(q));
}

template < class R >
inline
bool
y_equal(const Point_3<R> &p,
        const Point_3<R> &q)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return y_equal(static_cast<const RPoint_3&>(p),
	         static_cast<const RPoint_3&>(q));
}

template < class R >
inline
bool
z_equal(const Point_3<R> &p,
        const Point_3<R> &q)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return z_equal(static_cast<const RPoint_3&>(p),
	         static_cast<const RPoint_3&>(q));
}

template < class R >
inline
bool
coplanar(const Point_3<R> &p,
         const Point_3<R> &q,
         const Point_3<R> &r,
         const Point_3<R> &s)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return coplanar(static_cast<const RPoint_3&>(p),
	          static_cast<const RPoint_3&>(q),
                  static_cast<const RPoint_3&>(r),
		  static_cast<const RPoint_3&>(s));
}

template < class R >
inline
bool
collinear(const Point_3<R> &p,
          const Point_3<R> &q,
          const Point_3<R> &r)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return (collinear(static_cast<const RPoint_3&>(p),
                    static_cast<const RPoint_3&>(q),
                    static_cast<const RPoint_3&>(r)));
}

template < class R >
inline
bool
are_ordered_along_line(const Point_3<R> &p,
                       const Point_3<R> &q,
                       const Point_3<R> &r)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return (are_ordered_along_line(static_cast<const RPoint_3&>(p),
                                 static_cast<const RPoint_3&>(q),
                                 static_cast<const RPoint_3&>(r)));
}

template < class R >
inline
bool
collinear_are_ordered_along_line(const Point_3<R> &p,
                                 const Point_3<R> &q,
                                 const Point_3<R> &r)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return collinear_are_ordered_along_line(static_cast<const RPoint_3&>(p),
                                          static_cast<const RPoint_3&>(q),
                                          static_cast<const RPoint_3&>(r));
}

template < class R >
inline
bool
are_strictly_ordered_along_line(const Point_3<R> &p,
                                const Point_3<R> &q,
                                const Point_3<R> &r)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return (are_strictly_ordered_along_line(static_cast<const RPoint_3&>(p),
                                          static_cast<const RPoint_3&>(q),
                                          static_cast<const RPoint_3&>(r)));
}

template < class R >
inline
bool
collinear_are_strictly_ordered_along_line(const Point_3<R> &p,
                                          const Point_3<R> &q,
                                          const Point_3<R> &r)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return collinear_are_strictly_ordered_along_line(
	       static_cast<const RPoint_3&>(p),
	       static_cast<const RPoint_3&>(q),
	       static_cast<const RPoint_3&>(r));
}

template < class R >
inline
Orientation
orientation(const Point_3<R> &p,
            const Point_3<R> &q,
            const Point_3<R> &r,
            const Point_3<R> &s)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return orientation(static_cast<const RPoint_3&>(p),
	             static_cast<const RPoint_3&>(q),
                     static_cast<const RPoint_3&>(r),
		     static_cast<const RPoint_3&>(s));
}

template <class R >
inline
Bounded_side
side_of_bounded_sphere( const Point_3<R> &p,
                        const Point_3<R> &q,
                        const Point_3<R> &r,
                        const Point_3<R> &s,
                        const Point_3<R> &test)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return side_of_bounded_sphere(static_cast<const RPoint_3&>(p),
                                static_cast<const RPoint_3&>(q),
                                static_cast<const RPoint_3&>(r),
                                static_cast<const RPoint_3&>(s),
                                static_cast<const RPoint_3&>(test));
}

template <class R >
inline
Oriented_side
side_of_oriented_sphere( const Point_3<R> &p,
                         const Point_3<R> &q,
                         const Point_3<R> &r,
                         const Point_3<R> &s,
                         const Point_3<R> &test)
{
  typedef typename  R::Point_3_base  RPoint_3;
  return side_of_oriented_sphere(static_cast<const RPoint_3&>(p),
                                 static_cast<const RPoint_3&>(q),
                                 static_cast<const RPoint_3&>(r),
                                 static_cast<const RPoint_3&>(s),
                                 static_cast<const RPoint_3&>(test));
}

CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_ON_POINTS_3_H
