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

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#include <CGAL/TriangleH3.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/Triangle_3.h>
#endif

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

  bool                operator==(const CGAL::Triangle_3<R>& t) const
                      { return RTriangle_3::operator==(t); }
  bool                operator!=(const CGAL::Triangle_3<R>& t) const
                      { return !(*this == t); }
  CGAL::Plane_3<R>     supporting_plane() const
                      {
                        return
                        CGAL::Plane_3<R>(
                            RTriangle_3::supporting_plane());
                      }
  CGAL::Triangle_3<R>  transform(
                      const CGAL::Aff_transformation_3<R>& t) const
                      {
                        return
                        CGAL::Triangle_3<R>(RTriangle_3::transform( t ));
                      }
  bool                has_on(const CGAL::Point_3<R>& p) const
                      { return RTriangle_3::has_on(p); }
  bool                is_degenerate() const
                      { return RTriangle_3::is_degenerate(); }
  CGAL::Point_3<R>     vertex(int i) const
                      { return RTriangle_3::vertex(i); }
  CGAL::Point_3<R>     operator[](int i) const
                      { return vertex(i); }
  Bbox_3         bbox() const
                      {
                        return vertex(0).bbox()
                             + vertex(1).bbox()
                             + vertex(2).bbox();
                      }
};

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLE_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Triangle_3<R>& t)
{
  typedef typename  R::Triangle_3_base  RTriangle_3;
  return os << (const RTriangle_3& )t;
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLE_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_3
template < class R >
std::istream&
operator>>(std::istream& is, Triangle_3<R>& t)
{
  typedef typename  R::Triangle_3_base  RTriangle_3;
  return is >> (RTriangle_3& )t;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLE_3

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGLE_3_H
