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
// file          : Line_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_LINE_3_H
#define CGAL_LINE_3_H

#include <CGAL/Segment_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Plane_3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Line_3 : public R_::Line_3_base
{
public:
  typedef          R_                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Line_3_base  RLine_3;

  Line_3() : RLine_3()
  {}

  Line_3(const CGAL::Line_3<R>  & l)
      : RLine_3( static_cast<const RLine_3&>(l))
  {}

  Line_3(const CGAL::Point_3<R> & p, const CGAL::Point_3<R> & q)
      : RLine_3(p,q)
  {}

  // conversion impl -> interface class
  Line_3(const RLine_3&  l) : RLine_3(l)
  {}

  Line_3(const CGAL::Segment_3<R> & s) : RLine_3( s )
  {}

  Line_3(const CGAL::Ray_3<R> & r) : RLine_3( r )
  {}

  Line_3(const CGAL::Point_3<R> & p, const CGAL::Direction_3<R> & d)
      : RLine_3( p, d )
  {}

  CGAL::Plane_3<R>     perpendicular_plane(const CGAL::Point_3<R> & p) const
  { return RLine_3::perpendicular_plane(p); }

  CGAL::Line_3<R>      opposite() const
  { return RLine_3::opposite(); }

  CGAL::Point_3<R>     point() const
  { return RLine_3::point(); }

  CGAL::Point_3<R>     point(int i) const
  { return RLine_3::point(i); }

  CGAL::Point_3<R>     projection(const CGAL::Point_3<R>& p) const
  { return RLine_3::projection(p); }

  CGAL::Direction_3<R> direction() const
  { return RLine_3::direction(); }

  CGAL::Line_3<R>      transform(const CGAL::Aff_transformation_3<R> & t) const
  { return RLine_3::transform(t); }
};

#ifndef CGAL_NO_OSTREAM_INSERT_LINE_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Line_3<R>& l)
{
  typedef typename  R::Line_3_base  RLine_3;
  return os << static_cast<const RLine_3&>(l);
}
#endif // CGAL_NO_OSTREAM_INSERT_LINE_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINE_3
template < class R >
std::istream&
operator>>(std::istream & is, Line_3<R> & p)
{
  typedef typename  R::Line_3_base  RLine_3;
  is >> static_cast<RLine_3&>(p);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINE_3

CGAL_END_NAMESPACE

#endif // CGAL_LINE_3_H
