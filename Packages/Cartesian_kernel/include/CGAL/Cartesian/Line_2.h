// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/Line_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_LINE_2_H
#define CGAL_CARTESIAN_LINE_2_H

#include <CGAL/Threetuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class LineC2
  : public R_::template Handle<Threetuple<typename R_::FT> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Ray_2                Ray_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Threetuple<FT>	                   rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef R_                                     R;

  LineC2()
    : base(rep()) {}

  LineC2(const Point_2 &p, const Point_2 &q)
    : base(line_from_points(p, q)) {}

  LineC2(const FT &a, const FT &b, const FT &c)
    : base(rep(a, b, c)) {}

  LineC2(const Segment_2 &s)
    : base(line_from_points(s.source(), s.target())) {}

  LineC2(const Ray_2 &r)
    : base(line_from_points(r.point(0), r.point(1))) {}

  LineC2(const Point_2 &p, const Direction_2 &d)
    : base(line_from_point_direction(p, d)) {}

  bool            operator==(const LineC2 &l) const;
  bool            operator!=(const LineC2 &l) const;

  const FT & a() const
  {
      return Ptr()->e0;
  }
  const FT & b() const
  {
      return Ptr()->e1;
  }
  const FT & c() const
  {
      return Ptr()->e2;
  }

  FT              x_at_y(const FT &y) const;
  FT              y_at_x(const FT &x) const;

  Line_2          perpendicular(const Point_2 &p) const;
  Line_2          opposite() const;
  Point_2         point(int i) const;

  Point_2         point() const;
  Point_2         projection(const Point_2 &p) const;

  Direction_2     direction() const;

  Oriented_side   oriented_side(const Point_2 &p) const;
  bool            has_on_boundary(const Point_2 &p) const;
  bool            has_on_positive_side(const Point_2 &p) const;
  bool            has_on_negative_side(const Point_2 &p) const;
  bool            has_on(const Point_2 &p) const { return has_on_boundary(p); }

  bool            is_horizontal() const;
  bool            is_vertical() const;
  bool            is_degenerate() const;

  Line_2          transform(const Aff_transformation_2 &t) const
  {
    return LineC2<R>(t.transform(point(0)),
                               t.transform(direction()));
  }
};

template < class R >
CGAL_KERNEL_INLINE
bool
LineC2<R>::operator==(const LineC2<R> &l) const
{
  if (identical(l))
      return true;
  return equal_line(*this, l);
}

template < class R >
inline
bool
LineC2<R>::operator!=(const LineC2<R> &l) const
{
  return !(*this == l);
}

template < class R >
inline
bool
LineC2<R>::is_horizontal() const
{ // FIXME : predicate
  return CGAL_NTS is_zero(a());
}

template < class R >
inline
bool
LineC2<R>::is_vertical() const
{ // FIXME : predicate
  return CGAL_NTS is_zero(b());
}

template < class R >
CGAL_KERNEL_INLINE
typename LineC2<R>::FT
LineC2<R>::x_at_y(const typename LineC2<R>::FT &y) const
{
  CGAL_kernel_precondition_msg( ! is_horizontal(),
    "Line::x_at_y(FT y) is undefined for horizontal line");
  return line_x_at_y(*this, y);
}

template < class R >
CGAL_KERNEL_INLINE
typename LineC2<R>::FT
LineC2<R>::y_at_x(const typename LineC2<R>::FT &x) const
{
  CGAL_kernel_precondition_msg( ! is_vertical(),
    "Line::y_at_x(FT x) is undefined for vertical line");
  return line_y_at_x(*this, x);
}

template < class R >
inline
typename LineC2<R>::Line_2
LineC2<R>::
perpendicular(const typename LineC2<R>::Point_2 &p) const
{
  return perpendicular_through_point(*this, p);
}

template < class R >
inline
typename LineC2<R>::Line_2
LineC2<R>::opposite() const
{
  return LineC2<R>( -a(), -b(), -c() );
}

template < class R >
CGAL_KERNEL_INLINE
typename LineC2<R>::Point_2
LineC2<R>::point(int i) const
{
  return line_get_point(*this, i);
}

template < class R >
CGAL_KERNEL_INLINE
typename LineC2<R>::Point_2
LineC2<R>::point() const
{
  return line_get_point(*this, 0);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename LineC2<R>::Point_2
LineC2<R>::
projection(const typename LineC2<R>::Point_2 &p) const
{
  return line_project_point(*this, p);
}

template < class R >
inline
typename LineC2<R>::Direction_2
LineC2<R>::direction() const
{
  return Direction_2( b(), -a() );
}

template < class R >
CGAL_KERNEL_INLINE
Oriented_side
LineC2<R>::
oriented_side(const typename LineC2<R>::Point_2 &p) const
{
  return side_of_oriented_line(*this, p);
}

template < class R >
inline
bool
LineC2<R>::
has_on_boundary(const typename LineC2<R>::Point_2 &p) const
{
  return oriented_side(p) == ON_ORIENTED_BOUNDARY;
}

template < class R >
inline
bool
LineC2<R>::
has_on_positive_side(const typename LineC2<R>::Point_2 &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
bool
LineC2<R>::
has_on_negative_side(const typename LineC2<R>::Point_2 &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
LineC2<R>::is_degenerate() const
{
  return is_horizontal() && is_vertical();
}

#ifndef CGAL_NO_OSTREAM_INSERT_LINEC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const LineC2<R> &l)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << l.a() << ' ' << l.b() << ' ' << l.c();
    case IO::BINARY :
        write(os, l.a());
        write(os, l.b());
        write(os, l.c());
        return os;
    default:
        return os << "LineC2(" << l.a() 
		  << ", " << l.b() << ", " << l.c() <<')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_LINEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_LINEC2
template < class R >
std::istream &
operator>>(std::istream &is, LineC2<R> &l)
{
    typename R::FT a, b, c;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> a >> b >> c;
        break;
    case IO::BINARY :
        read(is, a);
        read(is, b);
        read(is, c);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
	l = LineC2<R>(a, b, c);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINEC2

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_LINE_2_H
