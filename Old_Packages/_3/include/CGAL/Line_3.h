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
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_LINE_3_H
#define CGAL_LINE_3_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Line_3 : public R_::Line_3_base
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_3               Point_3;
  typedef typename R_::Ray_3                 Ray_3;
  typedef typename R_::Segment_3             Segment_3;
  typedef typename R_::Direction_3           Direction_3;
  typedef typename R_::Line_3_base  RLine_3;
public:
  typedef          R_                       R;

  Line_3() : RLine_3()
  {}

  Line_3(const CGAL::Line_3<R>  & l)
      : RLine_3( static_cast<const RLine_3&>(l)) {}

  Line_3(const Point_3 & p, const Point_3 & q)
      : RLine_3(p,q) {}

  // conversion impl -> interface class
  Line_3(const RLine_3& l)
      : RLine_3(l) {}

  Line_3(const Segment_3 & s)
      : RLine_3( s ) {}

  Line_3(const Ray_3 & r)
      : RLine_3( r ) {}

  Line_3(const Point_3 & p, const Direction_3 & d)
      : RLine_3( p, d ) {}
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
