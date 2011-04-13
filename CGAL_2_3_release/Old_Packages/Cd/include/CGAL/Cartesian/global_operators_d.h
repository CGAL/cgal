// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/global_operators_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_GLOBAL_OPERATORS_D_C
#define CGAL_CARTESIAN_GLOBAL_OPERATORS_D_C

#include <CGAL/Cartesian/redefine_names_d.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
PointCd<R CGAL_CTAG>
operator+(const PointCd<R CGAL_CTAG> &p, const VectorCd<R CGAL_CTAG> &v)
{
  PointCd<R CGAL_CTAG> result(v.dimension());
  // Note: to be able to access the non-const result.begin(),
  // this function must be a friend of PointCd
  std::transform(p.begin(), p.end(), v.begin(), result.begin(),
                 std::plus<typename R::FT>());
  return result;
}

template < class R >
inline
PointCd<R CGAL_CTAG>
operator-(const PointCd<R CGAL_CTAG> &p, const VectorCd<R CGAL_CTAG> &v)
{
  PointCd<R CGAL_CTAG> result(v.dimension());
  // Note: to be able to access the non-const result.begin(),
  // this function must be a friend of PointCd
  std::transform(p.begin(), p.end(), v.begin(), result.begin(),
                 std::minus<typename R::FT>());
  return result;
}

template < class R >
inline
PointCd<R CGAL_CTAG>
operator+(const Origin &, const VectorCd<R CGAL_CTAG> &v)
{
  return PointCd<R CGAL_CTAG>(v);
}

template < class R >
inline
PointCd<R CGAL_CTAG>
operator-(const Origin &, const VectorCd<R CGAL_CTAG> &v)
{
  return PointCd<R CGAL_CTAG>(-v);
}

template < class R >
inline
VectorCd<R CGAL_CTAG>
operator-(const PointCd<R CGAL_CTAG> &p, const PointCd<R CGAL_CTAG> &q)
{
  VectorCd<R CGAL_CTAG> result(p.dimension());
  // Note: to be able to access the non-const result.begin(),
  // this function must be a friend of VectorCd
  std::transform(p.begin(), p.end(), q.begin(), result.begin(),
                 std::minus<typename R::FT>());
  return result;
}

template < class R >
inline
VectorCd<R CGAL_CTAG>
operator-(const PointCd<R CGAL_CTAG> &p, const Origin &)
{
  return VectorCd<R CGAL_CTAG>(p);
}

template < class R >
inline
VectorCd<R CGAL_CTAG>
operator-(const Origin &, const PointCd<R CGAL_CTAG> &p)
{
  return - VectorCd<R CGAL_CTAG>(p);
}

template < class R >
CGAL_KERNEL_INLINE
VectorCd<R CGAL_CTAG>
operator*(const typename R::FT &c, const VectorCd<R CGAL_CTAG> &w)
{
  return w * c;
}

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_GLOBAL_OPERATORS_D_C
