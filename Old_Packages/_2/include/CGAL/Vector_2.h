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
// file          : Vector_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_VECTOR_2_H
#define CGAL_VECTOR_2_H

CGAL_BEGIN_NAMESPACE

class Null_vector;

template <class R_>
class Vector_2 : public R_::Kernel_base::Vector_2
{
  typedef typename R_::RT             RT;
  typedef typename R_::Point_2        Point_2;
  typedef typename R_::Direction_2    Direction_2;
  typedef typename R_::Kernel_base::Vector_2  RVector_2;
public:
  typedef  R_                        R;

  Vector_2() {}

  Vector_2(const CGAL::Vector_2<R> &v)
      : RVector_2(static_cast<const RVector_2&>(v)) {}

  Vector_2(const Point_2& a, const Point_2& b)
      : RVector_2(a, b) {}

  Vector_2(const RVector_2& v) : RVector_2(v) {}

  Vector_2(const Null_vector &v) : RVector_2(v) {}

  Vector_2(const RT &x, const RT &y) : RVector_2(x,y) {}

  Vector_2(const RT &x, const RT &y, const RT &w) : RVector_2(x,y,w) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_VECTOR_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Vector_2<R> &v)
{
  typedef typename  R::Kernel_base::Vector_2  RVector_2;
  return os << static_cast<const RVector_2&>(v);
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTOR_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTOR_2
template < class R >
std::istream &
operator>>(std::istream &is, Vector_2<R> &p)
{
  typedef typename  R::Kernel_base::Vector_2  RVector_2;
  return is >> static_cast<RVector_2&>(p);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTOR_2

CGAL_END_NAMESPACE

#endif // CGAL_VECTOR_2_H
