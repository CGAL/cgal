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
// file          : include/CGAL/Cartesian/Line_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_LINE_2_C
#define CGAL_CARTESIAN_LINE_2_C

#include <CGAL/Cartesian/constructions_on_lines_2.h>
#include <CGAL/Cartesian/predicates_on_lines_2.h>

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
CGAL_KERNEL_INLINE
void
LineC2<R CGAL_CTAG>::new_rep(const typename LineC2<R CGAL_CTAG>::FT &a,
                             const typename LineC2<R CGAL_CTAG>::FT &b,
			     const typename LineC2<R CGAL_CTAG>::FT &c)
{
  new ( static_cast< void*>(ptr)) Threetuple<FT>(a, b, c);
}

template < class R >
CGAL_KERNEL_INLINE
void
LineC2<R CGAL_CTAG>::new_rep(const typename LineC2<R CGAL_CTAG>::Point_2 &p,
                             const typename LineC2<R CGAL_CTAG>::Point_2 &q)
{
  LineC2<R CGAL_CTAG> l = line_from_points(p,q);
  new_rep(l.a(), l.b(), l.c());
}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2()
{
  new ( static_cast< void*>(ptr)) Threetuple<FT>();
}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2(const LineC2<R CGAL_CTAG>  &l)
  : Handle_for<Threetuple<typename R::FT> >(l)
{}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2(const typename LineC2<R CGAL_CTAG>::Point_2 &p,
                            const typename LineC2<R CGAL_CTAG>::Point_2 &q)
{
  new_rep(p,q);
}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2(const typename LineC2<R CGAL_CTAG>::FT &a,
                            const typename LineC2<R CGAL_CTAG>::FT &b,
			    const typename LineC2<R CGAL_CTAG>::FT &c)
{
  new_rep(a, b, c);
}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2(const typename LineC2<R CGAL_CTAG>::Segment_2 &s)
{
  new_rep(s.source(), s.target());
}

template < class R >
inline
LineC2<R CGAL_CTAG>::LineC2(const typename LineC2<R CGAL_CTAG>::Ray_2 &r)
{
  new_rep(r.point(0), r.point(1));
}

template < class R >
CGAL_KERNEL_INLINE
LineC2<R CGAL_CTAG>::LineC2(const typename LineC2<R CGAL_CTAG>::Point_2 &p,
                            const typename LineC2<R CGAL_CTAG>::Direction_2 &d)
{
  LineC2<R CGAL_CTAG> l = line_from_point_direction(p,d);
  new_rep(l.a(), l.b(), l.c());
}

template < class R >
inline
LineC2<R CGAL_CTAG>::~LineC2()
{}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool LineC2<R CGAL_CTAG>::operator==(const LineC2<R CGAL_CTAG> &l) const
{
  if ( ptr == l.ptr ) return true;
  return equal_line(*this,l);
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
typename LineC2<R CGAL_CTAG>::FT
LineC2<R CGAL_CTAG>::a() const
{
  return ptr->e0;
}

template < class R >
inline
typename LineC2<R CGAL_CTAG>::FT
LineC2<R CGAL_CTAG>::b() const
{
  return ptr->e1;
}

template < class R >
inline
typename LineC2<R CGAL_CTAG>::FT
LineC2<R CGAL_CTAG>::c() const
{
  return ptr->e2;
}

template < class R >
inline
bool
LineC2<R CGAL_CTAG>::is_horizontal() const
{
  return a() == FT(0);
}

template < class R >
inline
bool
LineC2<R CGAL_CTAG>::is_vertical() const
{
  return b() == FT(0);
}

template < class R >
CGAL_KERNEL_INLINE
typename LineC2<R CGAL_CTAG>::FT
LineC2<R CGAL_CTAG>::x_at_y(const typename LineC2<R CGAL_CTAG>::FT &y) const
{
  CGAL_kernel_precondition_msg( a() != FT(0),
  "Line::x_at_y(FT y) is undefined for horizontal line");
  return line_x_at_y(*this,y);
}

template < class R >
CGAL_KERNEL_INLINE
typename LineC2<R CGAL_CTAG>::FT
LineC2<R CGAL_CTAG>::y_at_x(const typename LineC2<R CGAL_CTAG>::FT &x) const
{
  CGAL_kernel_precondition_msg( b() != FT(0),
  "Line::y_at_x(FT x) is undefined for vertical line");
  return line_y_at_x(*this,x);
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
{
  return LineC2<R CGAL_CTAG>( -a(), -b(), -c() );
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename LineC2<R CGAL_CTAG>::Point_2
LineC2<R CGAL_CTAG>::point(int i) const
{
  return line_get_point(*this,i);
}

template < class R >
CGAL_KERNEL_INLINE
typename LineC2<R CGAL_CTAG>::Point_2
LineC2<R CGAL_CTAG>::point() const
{
  return line_get_point(*this,0);
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
{
  return Direction_2( b(), -a() );
}

template < class R >
CGAL_KERNEL_INLINE
Oriented_side
LineC2<R CGAL_CTAG>::
oriented_side(const typename LineC2<R CGAL_CTAG>::Point_2 &p) const
{
  return side_of_oriented_line(*this,p);
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
  return a() == FT(0) && b() == FT(0);
}

template < class R >
inline
LineC2<R CGAL_CTAG>
LineC2<R CGAL_CTAG>::
transform(const typename LineC2<R CGAL_CTAG>::Aff_transformation_2 &t) const
{
  return LineC2<R CGAL_CTAG>(t.transform(point(0)),
                             t.transform(direction()));
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
operator>>(std::istream &is, LineC2<R CGAL_CTAG> &p)
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
    p = LineC2<R CGAL_CTAG>(a, b, c);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_LINEC2

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_LINE_2_C
