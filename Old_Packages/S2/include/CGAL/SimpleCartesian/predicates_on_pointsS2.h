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
// release_date  : 2000, August 16
//
// source        : webS2/S2.lw
// file          : include/CGAL/SimpleCartesian/predicates_on_pointsS2.h
// package       : S2 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.6
// revision_date : 27 Jun 2000
// author(s)     : Stefan Schirra
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================


#ifndef CGAL_PREDICATES_ON_POINTSS2_H
#define CGAL_PREDICATES_ON_POINTSS2_H

#include <CGAL/SimpleCartesian/simple_cartesian_classes.h>
#include <CGAL/SimpleCartesian/PointS2.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
inline
bool
x_equal(const PointS2<FT> &p, const PointS2<FT> &q)
{ return p.x() == q.x(); }

template < class FT >
inline
bool
y_equal(const PointS2<FT> &p, const PointS2<FT> &q)
{ return p.y() == q.y(); }

template < class FT >
inline
bool
equal_xy(const PointS2<FT> &p, const PointS2<FT> &q)
{ return (p.x() == q.x()) && (p.y() == q.y()); }

template < class FT >
inline
bool
less_x(const PointS2<FT> &p, const PointS2<FT> &q)
{ return (p.x() < q.x()); }

template < class FT >
inline
bool
less_y(const PointS2<FT> &p, const PointS2<FT> &q)
{ return (p.y() < q.y()); }

template < class FT >
inline
Comparison_result
compare_x(const PointS2<FT> &p, const PointS2<FT> &q)
{ return CGAL_NTS compare(p.x(), q.x()); }

template < class FT >
inline
Comparison_result
compare_y(const PointS2<FT> &p, const PointS2<FT> &q)
{ return CGAL_NTS compare(p.y(), q.y()); }

template < class FT >
inline
Comparison_result
compare_deltax_deltay(const PointS2<FT>& p,
                      const PointS2<FT>& q,
                      const PointS2<FT>& r,
                      const PointS2<FT>& s)
{
  return compare_deltax_deltayC2(p.x(), q.x(), r.y(), s.y());
}

template < class FT >
inline
Comparison_result
compare_lexicographically_xy(const PointS2<FT> &p,
                             const PointS2<FT> &q)
{
  return compare_lexicographically_xyC2(p.x(),p.y(),q.x(),q.y());
}

template < class FT >
inline
bool
lexicographically_xy_smaller_or_equal(const PointS2<FT> &p,
                                      const PointS2<FT> &q)
{
  return ( !( compare_lexicographically_xy(p,q) == LARGER ) );
}

template < class FT >
inline
bool
lexicographically_xy_smaller(const PointS2<FT> &p,
                             const PointS2<FT> &q)
{
  return compare_lexicographically_xy(p,q) == SMALLER ;
}

template < class FT >
inline
Comparison_result
compare_lexicographically_yx(const PointS2<FT> &p,
                             const PointS2<FT> &q)
{
  return compare_lexicographically_xyC2(p.y(),p.x(),q.y(),q.x());
}


template < class FT >
inline
bool
lexicographically_yx_smaller_or_equal(const PointS2<FT> &p,
                                      const PointS2<FT> &q)
{
    return  !( compare_lexicographically_yx(p,q) == LARGER ) ;
}

template < class FT >
inline
bool
lexicographically_yx_smaller(const PointS2<FT> &p,
                             const PointS2<FT> &q)
{
    return compare_lexicographically_yx(p,q) == SMALLER ;
}

template < class FT >
inline
Orientation
orientation(const PointS2<FT> &p,
            const PointS2<FT> &q,
            const PointS2<FT> &r)
{
    return orientationC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y());
}

template < class FT >
inline
bool
collinear(const PointS2<FT> &p,
               const PointS2<FT> &q,
               const PointS2<FT> &r)
{
  return (orientation(p,q,r) == COLLINEAR);
}



template < class FT >
inline
bool
collinear_are_ordered_along_line(const PointS2<FT> &p,
                                 const PointS2<FT> &q,
                                 const PointS2<FT> &r)
{
  return collinear_are_ordered_along_lineC2
             (p.x(),p.y(),q.x(),q.y(),r.x(),r.y());
}


template < class FT >
inline
bool
are_ordered_along_line(const PointS2<FT> &p,
                       const PointS2<FT> &q,
                       const PointS2<FT> &r)
{
  if (!collinear(p, q, r)) { return false; }
  return collinear_are_ordered_along_line(p, q, r);
}

template < class FT >
inline
bool
collinear_are_strictly_ordered_along_line(const PointS2<FT> &p,
                                          const PointS2<FT> &q,
                                          const PointS2<FT> &r)
{
  return collinear_are_strictly_ordered_along_lineC2
               (p.x(),p.y(),q.x(),q.y(),r.x(),r.y());
}


template < class FT >
inline
bool
are_strictly_ordered_along_line(const PointS2<FT> &p,
                                const PointS2<FT> &q,
                                const PointS2<FT> &r)
{
  if (!collinear(p, q, r)) { return false; }
  return collinear_are_strictly_ordered_along_line(p, q, r);
}

template < class FT >
inline
bool
leftturn(const PointS2<FT> &p,
         const PointS2<FT> &q,
         const PointS2<FT> &r)
{
  return (orientation(p,q,r) == LEFTTURN );
}

template < class FT >
inline
bool
leftturn(const Origin &o,
         const PointS2<FT> &q,
         const PointS2<FT> &r)
{
   return (orientationC2(FT(0), FT(0), q.x(), q.y(), r.x(), r.y())
           == LEFTTURN );
}

template < class FT >
inline
bool
rightturn(const PointS2<FT> &p,
          const PointS2<FT> &q,
          const PointS2<FT> &r)
{
   return (orientationC2(p.x(), p.y(), q.x(), q.y(), r.x(), r.y())
           == RIGHTTURN);
}

template < class FT >
inline
bool
rightturn(const Origin &o,
          const PointS2<FT> &q,
          const PointS2<FT> &r)
{
   return (orientationC2(FT(0), FT(0), q.x(), q.y(), r.x(), r.y())
           == RIGHTTURN);
}

template <class FT >
inline
Oriented_side
side_of_oriented_circle(const PointS2<FT> &p,
                        const PointS2<FT> &q,
                        const PointS2<FT> &r,
                        const PointS2<FT> &test)
{
  return side_of_oriented_circleC2
             (p.x(),p.y(),q.x(),q.y(),r.x(),r.y(),test.x(),test.y());
}


template <class FT >
inline
Bounded_side
side_of_bounded_circle(const PointS2<FT> &p,
                       const PointS2<FT> &q,
                       const PointS2<FT> &r,
                       const PointS2<FT> &test)
{
  return side_of_bounded_circleC2
             (p.x(),p.y(),q.x(),q.y(),r.x(),r.y(),test.x(),test.y());
}


CGAL_END_NAMESPACE

#endif  // CGAL_PREDICATES_ON_POINTSS2_H
