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
// file          : Vector_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_VECTOR_3_H
#define CGAL_VECTOR_3_H

CGAL_BEGIN_NAMESPACE

class Null_vector;

template <class R_>
class Vector_3 : public R_::Vector_3_base
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Direction_3           Direction_3;
  typedef typename R_::Vector_3_base         RVector_3;
public:
  typedef          R_                       R;

  Vector_3()
  {}

  Vector_3(const CGAL::Vector_3<R>& v)
    : RVector_3( static_cast<const RVector_3&>(v) ) {}

  Vector_3(const Point_3& a, const Point_3& b)
    : RVector_3(a, b) {}

  Vector_3(const RVector_3& v)
      : RVector_3(v) {}

  Vector_3(const Null_vector& v)
      : RVector_3(v) {}

  Vector_3(const RT& x, const RT& y, const RT& z)
    : RVector_3(x, y, z) {}

  Vector_3(const RT& x, const RT& y, const RT& z, const RT& w)
    : RVector_3(x, y, z, w) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_VECTOR_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Vector_3<R>& v)
{
  typedef typename  R::Vector_3_base  RVector_3;
  return os << static_cast<const RVector_3&>(v);
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTOR_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTOR_3
template < class R >
std::istream&
operator>>(std::istream& is, Vector_3<R>& p)
{
  typedef typename  R::Vector_3_base  RVector_3;
  return is >> static_cast<RVector_3&>(p);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTOR_3

CGAL_END_NAMESPACE

#endif // CGAL_VECTOR_3_H
