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
// file          : include/CGAL/Cartesian/Ray_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_RAY_3_H
#define CGAL_CARTESIAN_RAY_3_H

#include <CGAL/Cartesian/redefine_names_3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class RayC3 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Ray_handle_3
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Ray_handle_3              Ray_handle_3_;
  typedef typename Ray_handle_3_::element_type   Ray_ref_3;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef RayC3<R CGAL_CTAG>                    Self;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Direction_3               Direction_3;
  typedef typename R::Line_3                    Line_3;
  typedef typename R::Aff_transformation_3      Aff_transformation_3;
#else
  typedef RayC3<R>                              Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Direction_3_base          Direction_3;
  typedef typename R::Line_3_base               Line_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  RayC3()
    : Ray_handle_3_(Ray_ref_3()) {}

  RayC3(const Point_3 &sp, const Point_3 &secondp)
    : Ray_handle_3_(Ray_ref_3(sp, secondp)) {}

  RayC3(const Point_3 &sp, const Direction_3 &d)
    : Ray_handle_3_(Ray_ref_3(sp, sp + d.to_vector())) {}

  bool        operator==(const Self &r) const;
  bool        operator!=(const Self &r) const;

  Point_3     start() const;
  Point_3     source() const
  {
      return Ptr()->e0;
  }
  Point_3     second_point() const
  {
      return Ptr()->e1;
  }
  Point_3     point(int i) const;

  Direction_3 direction() const;
  Line_3      supporting_line() const;
  Self        opposite() const;

  Self        transform(const Aff_transformation_3 &t) const
  {
    return Self(t.transform(source()), t.transform(second_point()));
  }

  bool        is_degenerate() const;
  bool        has_on(const Point_3 &p) const;
  bool        collinear_has_on(const Point_3 &p) const;
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
inline
bool
RayC3<R CGAL_CTAG>::operator==(const RayC3<R CGAL_CTAG> &r) const
{
    if (identical(r))
	return true;
    return source() == r.source() && direction() == r.direction();
}

template < class R >
inline
bool
RayC3<R CGAL_CTAG>::operator!=(const RayC3<R CGAL_CTAG> &r) const
{
  return !(*this == r);
}

template < class R >
inline
typename RayC3<R CGAL_CTAG>::Point_3
RayC3<R CGAL_CTAG>::start() const
{
  return source();
}

template < class R >
CGAL_KERNEL_INLINE
typename RayC3<R CGAL_CTAG>::Point_3
RayC3<R CGAL_CTAG>::point(int i) const
{
  CGAL_kernel_precondition( i >= 0 );
  if (i == 0) return source();
  if (i == 1) return second_point();
  return source() + FT(i) * (second_point() - source());
}

template < class R >
inline
typename RayC3<R CGAL_CTAG>::Direction_3
RayC3<R CGAL_CTAG>::direction() const
{
  return Direction_3( second_point() - source() );
}

template < class R >
inline
typename RayC3<R CGAL_CTAG>::Line_3
RayC3<R CGAL_CTAG>::supporting_line() const
{
  return Line_3(*this);
}

template < class R >
inline
RayC3<R CGAL_CTAG>
RayC3<R CGAL_CTAG>::opposite() const
{
  return RayC3<R CGAL_CTAG>( source(), - direction() );
}

template < class R >
bool
RayC3<R CGAL_CTAG>::
has_on(const typename RayC3<R CGAL_CTAG>::Point_3 &p) const
{
  return (p == source()) ||
         ( collinear(source(), p, second_point())
           && ( Direction_3(p - source()) == direction() ));
}

template < class R >
inline
bool
RayC3<R CGAL_CTAG>::is_degenerate() const
{
  return source() == second_point();
}

template < class R >
inline
bool
RayC3<R CGAL_CTAG>::
collinear_has_on(const typename RayC3<R CGAL_CTAG>::Point_3 &p) const
{
  CGAL_kernel_exactness_precondition( collinear(source(), p, second_point()) );

  Comparison_result cx = compare_x(source(), second_point());
  if (cx != EQUAL)
    return cx != compare_x(p, source());

  Comparison_result cy = compare_y(source(), second_point());
  if (cy != EQUAL)
    return cy != compare_y(p, source());

  Comparison_result cz = compare_z(source(), second_point());
  if (cz != EQUAL)
    return cz != compare_z(p, source());

  return true; // p == source()
}

#ifndef CGAL_NO_OSTREAM_INSERT_RAYC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const RayC3<R CGAL_CTAG> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r.start() << ' ' << r.direction();
    case IO::BINARY :
        return os<< r.start() << r.direction();
    default:
        return os << "RayC3(" << r.start() <<  ", " << r.direction() << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_RAYC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_RAYC3
template < class R >
std::istream &
operator>>(std::istream &is, RayC3<R CGAL_CTAG> &r)
{
    typename RayC3<R CGAL_CTAG>::Point_3 p;
    typename RayC3<R CGAL_CTAG>::Direction_3 d;

    is >> p >> d;

    if (is)
	r = RayC3<R CGAL_CTAG>(p, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_RAYC3

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_RAY_3_H
