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
// author(s)     : Andreas Fabri, Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_TETRAHEDRON_3_H
#define CGAL_TETRAHEDRON_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Tetrahedron_3 : public R_::Kernel_base::Tetrahedron_3
{
  typedef typename R_::Point_3             Point_3;
  typedef typename R_::Kernel_base::Tetrahedron_3  RTetrahedron_3;
public:
  typedef          R_                       R;

  Tetrahedron_3()
      : RTetrahedron_3() {}

  Tetrahedron_3(const CGAL::Tetrahedron_3<R>& t)
      : RTetrahedron_3(t) {}

  Tetrahedron_3(const RTetrahedron_3& t)
      : RTetrahedron_3(t) {}

  Tetrahedron_3(const Point_3& p,
                const Point_3& q,
                const Point_3& r,
                const Point_3& s)
    : RTetrahedron_3(p,q,r,s) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_TETRAHEDRON_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Tetrahedron_3<R>& t)
{
  typedef typename  R::Kernel_base::Tetrahedron_3  RTetrahedron_3;
  return os << static_cast<const RTetrahedron_3&>(t);
}
#endif // CGAL_NO_OSTREAM_INSERT_TETRAHEDRON_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRON_3
template < class R >
std::istream&
operator>>(std::istream& is, Tetrahedron_3<R>& t)
{
  typedef typename  R::Kernel_base::Tetrahedron_3  RTetrahedron_3;
  return is >> static_cast<RTetrahedron_3&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRON_3

CGAL_END_NAMESPACE

#endif  // CGAL_TETRAHEDRON_3_H
