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

CGAL_BEGIN_NAMESPACE

template <class R_>
class Line_2 : public R_::Line_2_base
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Segment_2             Segment_2;
  typedef typename R_::Ray_2                 Ray_2;
  typedef typename R_::Direction_2           Direction_2;
  typedef typename R_::Line_2_base  RLine_2;
public:
  typedef  R_   R;

  Line_2()
    : RLine_2() {}

  Line_2(const CGAL::Line_2<R>  &l)
    : RLine_2(static_cast<const RLine_2&>(l)) {}

  Line_2(const Point_2 &p, const Point_2 &q)
    : RLine_2(p,q) {}

  Line_2(const RT &a, const RT &b, const RT &c)
    : RLine_2(a,b,c) {}

  Line_2(const RLine_2& l)  // conversion impl -> interface class
    : RLine_2(l) {}

  Line_2(const Segment_2& s)
    : RLine_2(s) {}

  Line_2(const Ray_2& r)
    : RLine_2(r) {}

  Line_2(const Point_2 &p, const Direction_2 &d)
    : RLine_2(p,d) {}
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
