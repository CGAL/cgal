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

#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/Cartesian/predicates_on_lines_2.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class LineC2 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Line_handle_2
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Line_handle_2             Line_handle_2_;
  typedef typename Line_handle_2_::element_type  Line_ref_2;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef LineC2<R,Cartesian_tag>               Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Vector_2                  Vector_2;
  typedef typename R::Direction_2               Direction_2;
  typedef typename R::Ray_2                     Ray_2;
  typedef typename R::Triangle_2                Triangle_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Iso_rectangle_2           Iso_rectangle_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
  typedef typename R::Circle_2                  Circle_2;
#else
  typedef LineC2<R>                             Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Vector_2_base             Vector_2;
  typedef typename R::Direction_2_base          Direction_2;
  typedef typename R::Ray_2_base                Ray_2;
  typedef typename R::Triangle_2_base           Triangle_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Iso_rectangle_2_base      Iso_rectangle_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
  typedef typename R::Circle_2_base             Circle_2;
#endif

  LineC2()
    : Line_handle_2_(Line_ref_2()) {}

  LineC2(const Point_2 &p, const Point_2 &q)
    : Line_handle_2_(line_from_points(p, q)) {}

  LineC2(const FT &a, const FT &b, const FT &c)
    : Line_handle_2_(Line_ref_2(a, b, c)) {}

  LineC2(const Segment_2 &s)
    : Line_handle_2_(line_from_points(s.source(), s.target())) {}

  LineC2(const Ray_2 &r)
    : Line_handle_2_(line_from_points(r.point(0), r.point(1))) {}

  LineC2(const Point_2 &p, const Direction_2 &d)
    : Line_handle_2_(line_from_point_direction(p, d)) {}

  bool            operator==(const Self &l) const;
  bool            operator!=(const Self &l) const;

  FT a() const
  {
      return Ptr()->e0;
  }
  FT b() const
  {
      return Ptr()->e1;
  }
  FT c() const
  {
      return Ptr()->e2;
  }

  FT              x_at_y(const FT &y) const;
  FT              y_at_x(const FT &x) const;

  Self            perpendicular(const Point_2 &p) const;
  Self            opposite() const;
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

  Self            transform(const Aff_transformation_2 &t) const
  {
    return LineC2<R CGAL_CTAG>(t.transform(point(0)),
                               t.transform(direction()));
  }
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_INLINE
bool
LineC2<R CGAL_CTAG>::operator==(const LineC2<R CGAL_CTAG> &l) const
{
  if (identical(l))
      return true;
  return equal_line(*this, l);
}

template < class R >
inline
bool
LineC2<R CGAL_CTAG>::operator!=(const LineC2<R CGAL_CTAG> &l) const
{
  return !(*this == l);
}

template < class R >
inline
bool
LineC2<R CGAL_CTAG>::is_horizontal() const
{ // FIXME : predicate
  return CGAL_NTS is_zero(a());
}

template < class R >
inline
bool
LineC2<R CGAL_CTAG>::is_vertical() const
{ // FIXME : predicate
  return CGAL_NTS is_zero(b());
}

template < class R >
CGAL_KERNEL_INLINE
typename LineC2<R CGAL_CTAG>::FT
LineC2<R CGAL_CTAG>::x_at_y(const typename LineC2<R CGAL_CTAG>::FT &y) const
{
  CGAL_kernel_precondition_msg( ! is_horizontal(),
    "Line::x_at_y(FT y) is undefined for horizontal line");
  return line_x_at_y(*this, y);
}

template < class R >
CGAL_KERNEL_INLINE
typename LineC2<R CGAL_CTAG>::FT
LineC2<R CGAL_CTAG>::y_at_x(const typename LineC2<R CGAL_CTAG>::FT &x) const
{
  CGAL_kernel_precondition_msg( ! is_vertical(),
    "Line::y_at_x(FT x) is undefined for vertical line");
  return line_y_at_x(*this, x);
}

template < class R >
inline
LineC2<R CGAL_CTAG>
LineC2<R CGAL_CTAG>::
perpendicular(const typename LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return perpendicular_through_point(*this, p);
}

template < class R >
inline
LineC2<R CGAL_CTAG>
LineC2<R CGAL_CTAG>::opposite() const
{ // FIXME : construction
  return LineC2<R CGAL_CTAG>( -a(), -b(), -c() );
}

template < class R >
CGAL_KERNEL_INLINE
typename LineC2<R CGAL_CTAG>::Point_2
LineC2<R CGAL_CTAG>::point(int i) const
{
  return line_get_point(*this, i);
}

template < class R >
CGAL_KERNEL_INLINE
typename LineC2<R CGAL_CTAG>::Point_2
LineC2<R CGAL_CTAG>::point() const
{
  return line_get_point(*this, 0);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename LineC2<R CGAL_CTAG>::Point_2
LineC2<R CGAL_CTAG>::
projection(const typename LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return line_project_point(*this, p);
}

template < class R >
inline
typename LineC2<R CGAL_CTAG>::Direction_2
LineC2<R CGAL_CTAG>::direction() const
{ // FIXME : construction
  return Direction_2( b(), -a() );
}

template < class R >
CGAL_KERNEL_INLINE
Oriented_side
LineC2<R CGAL_CTAG>::
oriented_side(const typename LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return side_of_oriented_line(*this, p);
}

template < class R >
inline
bool
LineC2<R CGAL_CTAG>::
has_on_boundary(const typename LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return oriented_side(p) == ON_ORIENTED_BOUNDARY;
}

template < class R >
inline
bool
LineC2<R CGAL_CTAG>::
has_on_positive_side(const typename LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
bool
LineC2<R CGAL_CTAG>::
has_on_negative_side(const typename LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
LineC2<R CGAL_CTAG>::is_degenerate() const
{
  return is_horizontal() && is_vertical();
}

#ifndef CGAL_NO_OSTREAM_INSERT_LINEC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const LineC2<R CGAL_CTAG> &l)
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
operator>>(std::istream &is, LineC2<R CGAL_CTAG> &l)
{
    typename LineC2<R CGAL_CTAG>::FT a, b, c;
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
	l = LineC2<R CGAL_CTAG>(a, b, c);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINEC2

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_LINE_2_H
