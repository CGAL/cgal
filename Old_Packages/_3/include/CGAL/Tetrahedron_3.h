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

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/TetrahedronH3.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/Tetrahedron_3.h>
#endif

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
  bool               operator==(const CGAL::Tetrahedron_3<R>& t) const
                     { return RTetrahedron_3::operator==(t); }
  bool               operator!=(const CGAL::Tetrahedron_3<R>& t) const
                     { return !(*this == t); }
  Bbox_3             bbox() const
                     {
                       return vertex(0).bbox() + vertex(1).bbox()
                            + vertex(2).bbox() + vertex(3).bbox();
                     }
  FT                 volume() const
                     {  return RTetrahedron_3::volume(); }
  CGAL::Tetrahedron_3<R>
                     transform(const CGAL::Aff_transformation_3<R>& t) const
                     {
                       return
                       CGAL::Tetrahedron_3<R>(RTetrahedron_3::transform(t));
                     }
  Orientation   orientation() const
                     { return RTetrahedron_3::orientation(); }
  Oriented_side oriented_side(const CGAL::Point_3<R>& p) const
                     { return RTetrahedron_3::oriented_side(p); }
  bool               has_on_positive_side(const CGAL::Point_3<R>& p) const
                     { return oriented_side(p) == ON_POSITIVE_SIDE; }
  bool               has_on_negative_side(const CGAL::Point_3<R>& p) const
                     { return oriented_side(p) == ON_NEGATIVE_SIDE; }
  Bounded_side  bounded_side(const CGAL::Point_3<R>& p) const
                     { return RTetrahedron_3::bounded_side(p); }
  bool               has_on_boundary(const CGAL::Point_3<R>& p) const
                     { return bounded_side(p) == ON_BOUNDARY; }
  bool               has_on_bounded_side(const CGAL::Point_3<R>& p) const
                     { return bounded_side(p) == ON_BOUNDED_SIDE; }
  bool               has_on_unbounded_side(const CGAL::Point_3<R>& p) const
                     { return bounded_side(p) == ON_UNBOUNDED_SIDE; }
  bool               is_degenerate() const
                     { return RTetrahedron_3::is_degenerate(); }
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
