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
// release_date  : 2000, October 15
//
// source        : webS3/S3.lw
// file          : include/CGAL/SimpleCartesian/predicates_on_pointsS3.h
// package       : S3 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.7
// revision_date : 15 Oct 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@@mpi-sb.mpg.de>
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================

#ifndef CGAL_PREDICATES_ON_POINTSS3_H
#define CGAL_PREDICATES_ON_POINTSS3_H

#include <CGAL/predicates/kernel_ftC3.h>

CGAL_BEGIN_NAMESPACE

template < class FT >
Comparison_result
compare_lexicographically_xyz(const PointS3<FT>& p, const PointS3<FT>& q)
{
  return compare_lexicographically_xyzC3(p.x(),p.y(),p.z(),
                                         q.x(),q.y(),q.z());
}

template < class FT >
Comparison_result
compare_lexicographically_xy(const PointS3<FT>& p, const PointS3<FT>& q)
{
  return compare_lexicographically_xyC2(p.x(),p.y(),
                                        q.x(),q.y());
}

template < class FT >
bool
lexicographically_xyz_smaller_or_equal(const PointS3<FT>& p,
                                       const PointS3<FT>& q)
{ return compare_lexicographically_xyz(p,q) != LARGER; }

template < class FT >
bool
lexicographically_xyz_smaller(const PointS3<FT>& p, const PointS3<FT>& q)
{ return compare_lexicographically_xyz(p,q) == SMALLER; }

template < class FT >
inline
bool
lexicographically_xy_smaller(const PointS3<FT>& p, const PointS3<FT>& q)
{ return compare_lexicographically_xy(p,q) == SMALLER; }

template < class FT >
inline
bool
x_equal(const PointS3<FT>& p, const PointS3<FT>& q)
{ return p.x() == q.x(); }


template < class FT >
inline
bool
y_equal(const PointS3<FT>& p, const PointS3<FT>& q)
{ return p.y() == q.y(); }

template < class FT >
inline
bool
z_equal(const PointS3<FT>& p, const PointS3<FT>& q)
{ return p.z() == q.z(); }

template < class FT >
inline
bool
equal_xy(const PointS3<FT>& p, const PointS3<FT>& q)
{ return p.x() == q.x() && p.y() == q.y(); }

template < class FT >
inline
bool
equal_xyz(const PointS3<FT>& p, const PointS3<FT>& q)
{ return p.x() == q.x() && p.y() == q.y() && p.z() == q.z(); }

template < class FT >
inline
Comparison_result
compare_x(const PointS3<FT>& p, const PointS3<FT>& q)
{ return CGAL_NTS compare(p.x(), q.x()); }


template < class FT >
inline
Comparison_result
compare_y(const PointS3<FT>& p, const PointS3<FT>& q)
{ return CGAL_NTS compare(p.y(), q.y()); }


template < class FT >
inline
Comparison_result
compare_z(const PointS3<FT>& p, const PointS3<FT>& q)
{ return CGAL_NTS compare(p.z(), q.z()); }

template < class FT >
inline
bool
less_x(const PointS3<FT>& p, const PointS3<FT>& q)
{ return p.x() < q.x(); }

template < class FT >
inline
bool
less_y(const PointS3<FT>& p, const PointS3<FT>& q)
{ return p.y() < q.y(); }

template < class FT >
inline
bool
less_z(const PointS3<FT>& p, const PointS3<FT>& q)
{ return p.z() < q.z(); }


template < class FT >
bool collinear(const PointS3<FT>& p,
               const PointS3<FT>& q,
               const PointS3<FT>& r)
{
  return collinearC3(p.x(), p.y(), p.z(),
                     q.x(), q.y(), q.z(),
                     r.x(), r.y(), r.z());
}

template < class FT >
inline
Orientation
orientation(const PointS3<FT>& p,
            const PointS3<FT>& q,
            const PointS3<FT>& r,
            const PointS3<FT>& s)
{
  return orientationC3(p.x(), p.y(), p.z(),
                       q.x(), q.y(), q.z(),
                       r.x(), r.y(), r.z(),
                       s.x(), s.y(), s.z());
}

template < class FT >
inline
bool
coplanar(const PointS3<FT>& p,
         const PointS3<FT>& q,
         const PointS3<FT>& r,
         const PointS3<FT>& s)
{
  return orientation(p, q, r, s) == COLLINEAR;
}

template < class FT>
inline
bool
are_positive_oriented(const PointS3<FT>& p,
                      const PointS3<FT>& q,
                      const PointS3<FT>& r,
                      const PointS3<FT>& s)
{
  return orientation(p,q,r,s) == POSITIVE;
}

template < class FT>
inline
bool
are_negative_oriented(const PointS3<FT>& p,
                      const PointS3<FT>& q,
                      const PointS3<FT>& r,
                      const PointS3<FT>& s)
{
  return orientation(p,q,r,s) == NEGATIVE;
}

template < class FT >
inline
bool
are_ordered_along_line(const PointS3<FT>& p,
                       const PointS3<FT>& q,
                       const PointS3<FT>& r)
{
  return (collinear(p, q, r)) ? collinear_are_ordered_along_line(p, q, r)
                              : false;
}

template < class FT >
inline
bool
collinear_are_ordered_along_line(const PointS3<FT>& p,
                                 const PointS3<FT>& q,
                                 const PointS3<FT>& r)
{
  CGAL_kernel_exactness_precondition( collinear(p, q, r) );
  return collinear_are_ordered_along_lineC3(p.x(),p.y(),p.z(),
                                            q.x(),q.y(),q.z(),
                                            r.x(),r.y(),r.z());
}


template < class FT >
inline
bool
are_strictly_ordered_along_line(const PointS3<FT>& p,
                                const PointS3<FT>& q,
                                const PointS3<FT>& r)
{
  return (collinear(p, q, r))
         ? collinear_are_strictly_ordered_along_line(p, q, r)
         : false;
}


template < class FT >
inline
bool
collinear_are_strictly_ordered_along_line(const PointS3<FT>& p,
                                          const PointS3<FT>& q,
                                          const PointS3<FT>& r)
{
  CGAL_kernel_exactness_precondition( collinear(p, q, r) );
  return collinear_are_strictly_ordered_along_lineC3(p.x(),p.y(),p.z(),
                                                     q.x(),q.y(),q.z(),
                                                     r.x(),r.y(),r.z());
}


template <class FT>
Orientation
coplanar_orientation(const PointS3<FT>& p,
                     const PointS3<FT>& q,
                     const PointS3<FT>& r,
                     const PointS3<FT>& s)
{
  CGAL_kernel_exactness_precondition( ! collinear(p, q, r) );
  CGAL_kernel_exactness_precondition( coplanar(s, p, q, r) );
  return coplanar_orientationC3(p.x(), p.y(), p.z(),
                                q.x(), q.y(), q.z(),
                                r.x(), r.y(), r.z(),
                                s.x(), s.y(), s.z());
}


template <class FT >
Oriented_side
side_of_oriented_sphere(const PointS3<FT>& p,
                        const PointS3<FT>& q,
                        const PointS3<FT>& r,
                        const PointS3<FT>& s,
                        const PointS3<FT>& test)
{
  return side_of_oriented_sphereC3(p.x(),p.y(),p.z(),
                                   q.x(),q.y(),q.z(),
                                   r.x(),r.y(),r.z(),
                                   s.x(),s.y(),s.z(),
                                   test.x(),test.y(),test.z());
}


template <class FT >
Bounded_side
side_of_bounded_sphere(const PointS3<FT>& p,
                       const PointS3<FT>& q,
                       const PointS3<FT>& r,
                       const PointS3<FT>& s,
                       const PointS3<FT>& test)
{
  return side_of_bounded_sphereC3(p.x(),p.y(),p.z(),
                                  q.x(),q.y(),q.z(),
                                  r.x(),r.y(),r.z(),
                                  s.x(),s.y(),s.z(),
                                  test.x(),test.y(),test.z());
}


CGAL_END_NAMESPACE

#endif // CGAL_PREDICATES_ON_POINTSS3_H
