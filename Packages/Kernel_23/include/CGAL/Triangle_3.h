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
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_TRIANGLE_3_H
#define CGAL_TRIANGLE_3_H

#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Triangle_3 : public R_::Triangle_3_base
{
public:
  typedef          R_                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Triangle_3_base  RTriangle_3;

  Triangle_3() : RTriangle_3()
  {}

  Triangle_3(const CGAL::Triangle_3<R>& t) : RTriangle_3(t)
  {}

  Triangle_3(const RTriangle_3&  t) : RTriangle_3(t)
  {}

  Triangle_3(const CGAL::Point_3<R>& p,
                  const CGAL::Point_3<R>& q,
                  const CGAL::Point_3<R>& r)
    : RTriangle_3(p,q,r)
  {}

  CGAL::Plane_3<R>     supporting_plane() const
  {
      return CGAL::Plane_3<R>(RTriangle_3::supporting_plane());
  }

  CGAL::Triangle_3<R>  transform(const CGAL::Aff_transformation_3<R>& t) const
  {
      return CGAL::Triangle_3<R>(RTriangle_3::transform( t ));
  }

  CGAL::Point_3<R>     vertex(int i) const
  { return RTriangle_3::vertex(i); }

  CGAL::Point_3<R>     operator[](int i) const
  { return vertex(i); }
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
