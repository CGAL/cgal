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
// file          : Triangle_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_TRIANGLE_3_H
#define CGAL_TRIANGLE_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Triangle_3 : public R_::Triangle_3_base
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Triangle_3_base  RTriangle_3;
public:
  typedef          R_                       R;

  Triangle_3()
      : RTriangle_3() {}

  Triangle_3(const CGAL::Triangle_3<R>& t)
      : RTriangle_3(t) {}

  Triangle_3(const RTriangle_3& t)
      : RTriangle_3(t) {}

  Triangle_3(const Point_3& p, const Point_3& q, const Point_3& r)
    : RTriangle_3(p,q,r) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLE_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Triangle_3<R>& t)
{
  typedef typename  R::Triangle_3_base  RTriangle_3;
  return os << static_cast<const RTriangle_3&>(t);
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLE_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_3
template < class R >
std::istream&
operator>>(std::istream& is, Triangle_3<R>& t)
{
  typedef typename  R::Triangle_3_base  RTriangle_3;
  return is >> static_cast<RTriangle_3&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_3

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGLE_3_H
