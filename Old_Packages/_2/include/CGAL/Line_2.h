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
// file          : Line_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================

#ifndef CGAL_LINE_2_H
#define CGAL_LINE_2_H

#include <CGAL/Segment_2.h>
#include <CGAL/Ray_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Line_2 : public R_::Line_2_base
{
public:
  typedef  R_   R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Line_2_base  RLine_2;

  Line_2()
    : RLine_2()
  {}

  Line_2(const CGAL::Line_2<R>  &l)
    : RLine_2(static_cast<const RLine_2&>(l))
  {}

  Line_2(const CGAL::Point_2<R> &p, const CGAL::Point_2<R> &q)
    : RLine_2(p,q)
  {}

  Line_2(const RT &a, const RT &b, const RT &c)
    : RLine_2(a,b,c)
  {}

  Line_2(const RLine_2& l)  // conversion impl -> interface class
    : RLine_2(l)
  {}

  Line_2(const CGAL::Segment_2<R>& s)
    : RLine_2(s)
  {}

  Line_2(const CGAL::Ray_2<R>& r)
    : RLine_2(r)
  {}

  Line_2(const CGAL::Point_2<R> &p, const CGAL::Direction_2<R> &d)
    : RLine_2(p,d)
  {}

  CGAL::Line_2<R> perpendicular(const CGAL::Point_2<R> &p) const
  {
    return RLine_2::perpendicular(p);
  }

  CGAL::Line_2<R> opposite() const
  {
    return RLine_2::opposite();
  }

  CGAL::Point_2<R> point(int i) const
  {
    return RLine_2::point(i);
  }

  CGAL::Point_2<R> projection(const CGAL::Point_2<R> &p) const
  {
    return RLine_2::projection(p);
  }

  CGAL::Point_2<R> point() const
  {
    return RLine_2::point();
  }

  CGAL::Direction_2<R> direction() const
  {
    return RLine_2::direction();
  }

  CGAL::Line_2<R> transform(const CGAL::Aff_transformation_2<R> &t) const
  {
    return  RLine_2::transform(t);
  }
};

#ifndef CGAL_NO_OSTREAM_INSERT_LINE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Line_2<R> &l)
{
  typedef typename  R::Line_2_base  RLine_2;
  return os << static_cast<const RLine_2&>(l);
}
#endif // CGAL_NO_OSTREAM_INSERT_LINE_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINE_2
template < class R >
std::istream &
operator>>(std::istream &is, Line_2<R> &p)
{
  typedef typename  R::Line_2_base  RLine_2;
  return is >> static_cast<RLine_2&>(p);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINE_2

CGAL_END_NAMESPACE

#endif  // CGAL_LINE_2_H
