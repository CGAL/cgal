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
// file          : Tetrahedron_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_TETRAHEDRON_3_H
#define CGAL_TETRAHEDRON_3_H

#include <CGAL/Plane_3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Tetrahedron_3 : public R_::Tetrahedron_3_base
{
public:
  typedef          R_                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Tetrahedron_3_base  RTetrahedron_3;

  Tetrahedron_3() : RTetrahedron_3()
  {}
  Tetrahedron_3(const CGAL::Tetrahedron_3<R>& t) : RTetrahedron_3(t)
  {}
  Tetrahedron_3(const RTetrahedron_3&  t) : RTetrahedron_3(t)
  {}
  Tetrahedron_3(const CGAL::Point_3<R>& p,
                     const CGAL::Point_3<R>& q,
                     const CGAL::Point_3<R>& r,
                     const CGAL::Point_3<R>& s)
    : RTetrahedron_3(p,q,r,s)
  {}

  CGAL::Point_3<R>    vertex(int i) const
                     { return RTetrahedron_3::vertex(i); }
  CGAL::Point_3<R>    operator[](int i) const
                     { return vertex(i); }
  CGAL::Tetrahedron_3<R>
                     transform(const CGAL::Aff_transformation_3<R>& t) const
                     {
                       return
                       CGAL::Tetrahedron_3<R>(RTetrahedron_3::transform(t));
                     }
};

#ifndef CGAL_NO_OSTREAM_INSERT_TETRAHEDRON_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Tetrahedron_3<R>& t)
{
  typedef typename  R::Tetrahedron_3_base  RTetrahedron_3;
  return os << static_cast<const RTetrahedron_3&>(t);
}
#endif // CGAL_NO_OSTREAM_INSERT_TETRAHEDRON_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRON_3
template < class R >
std::istream&
operator>>(std::istream& is, Tetrahedron_3<R>& t)
{
  typedef typename  R::Tetrahedron_3_base  RTetrahedron_3;
  return is >> static_cast<RTetrahedron_3&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRON_3

CGAL_END_NAMESPACE

#endif  // CGAL_TETRAHEDRON_3_H
