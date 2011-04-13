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
// file          : include/CGAL/Cartesian/Plane_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_PLANE_3_H
#define CGAL_CARTESIAN_PLANE_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Cartesian/predicates_on_planes_3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class PlaneC3 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Plane_handle_3
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Plane_handle_3            Plane_handle_3_;
  typedef typename Plane_handle_3_::element_type Plane_ref_3;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef PlaneC3<R,Cartesian_tag>              Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Vector_3                  Vector_3;
  typedef typename R::Direction_3               Direction_3;
  typedef typename R::Line_3                    Line_3;
  typedef typename R::Ray_3                     Ray_3;
  typedef typename R::Segment_3                 Segment_3;
  typedef typename R::Aff_transformation_3      Aff_transformation_3;
#else
  typedef PlaneC3<R>                            Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Vector_3_base             Vector_3;
  typedef typename R::Direction_3_base          Direction_3;
  typedef typename R::Line_3_base               Line_3;
  typedef typename R::Ray_3_base                Ray_3;
  typedef typename R::Segment_3_base            Segment_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  PlaneC3()
    : Plane_handle_3_(Plane_ref_3()) {}

  PlaneC3(const Point_3 &p, const Point_3 &q, const Point_3 &r)
    : Plane_handle_3_(plane_from_points(p, q, r)) {}

  PlaneC3(const Point_3 &p, const Direction_3 &d)
    : Plane_handle_3_(plane_from_point_direction(p, d)) {}

  PlaneC3(const Point_3 &p, const Vector_3 &v)
    : Plane_handle_3_(plane_from_point_direction(p, v.direction())) {}

  PlaneC3(const FT &a, const FT &b, const FT &c, const FT &d)
    : Plane_handle_3_(Plane_ref_3(a, b, c, d)) {}

  PlaneC3(const Line_3 &l, const Point_3 &p)
    : Plane_handle_3_(plane_from_points(l.point(),
	                               l.point()+l.direction().to_vector(),
				       p)) {}

  PlaneC3(const Segment_3 &s, const Point_3 &p)
    : Plane_handle_3_(plane_from_points(s.start(), s.end(), p)) {}

  PlaneC3(const Ray_3 &r, const Point_3 &p)
    : Plane_handle_3_(plane_from_points(r.start(), r.second_point(), p)) {}

  bool         operator==(const Self &p) const;
  bool         operator!=(const Self &p) const;

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
  FT d() const
  {
      return Ptr()->e3;
  }

  Line_3       perpendicular_line(const Point_3 &p) const;
  Self         opposite() const;

  Point_3      point() const;
  Point_3      projection(const Point_3 &p) const;
  Vector_3     orthogonal_vector() const;
  Direction_3  orthogonal_direction() const;
  Vector_3     base1() const;
  Vector_3     base2() const;

  Point_3      to_plane_basis(const Point_3 &p) const;

  Point_2      to_2d(const Point_3 &p) const;
  Point_3      to_3d(const Point_2 &p) const;

  Self         transform(const Aff_transformation_3 &t) const
  {
    if (t.is_even())
      return Self(t.transform(point()),
                 t.transpose().inverse().transform(orthogonal_direction()));
    else
      return Self( t.transform(point()),
               - t.transpose().inverse().transform(orthogonal_direction()));
  }

  Oriented_side oriented_side(const Point_3 &p) const;
#ifndef CGAL_NO_DEPRECATED_CODE
  bool         has_on_boundary(const Point_3 &p) const;
  bool         has_on_boundary(const Line_3 &l) const;
#endif // CGAL_NO_DEPRECATED_CODE
  bool         has_on_positive_side(const Point_3 &l) const;
  bool         has_on_negative_side(const Point_3 &l) const;
  bool         has_on(const Point_3 &p) const;
  bool         has_on(const Line_3 &l) const;

  bool         is_degenerate() const;
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_INLINE
bool
PlaneC3<R CGAL_CTAG>::operator==(const PlaneC3<R CGAL_CTAG> &p) const
{
  if (identical(p))
      return true;
  return equal_plane(*this, p);
}

template < class R >
inline
bool
PlaneC3<R CGAL_CTAG>::operator!=(const PlaneC3<R CGAL_CTAG> &p) const
{
  return !(*this == p);
}

template < class R >
inline
typename PlaneC3<R CGAL_CTAG>::Point_3
PlaneC3<R CGAL_CTAG>::point() const
{
  return point_on_plane(*this);
}

template < class R >
inline
typename PlaneC3<R CGAL_CTAG>::Point_3
PlaneC3<R CGAL_CTAG>::
projection(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return CGAL::projection_plane(p, *this); // FIXME : CGAL:: needed ?
}

template < class R >
inline
typename PlaneC3<R CGAL_CTAG>::Vector_3
PlaneC3<R CGAL_CTAG>::orthogonal_vector() const
{
  return Vector_3(a(), b(), c());
}

template < class R >
inline
typename PlaneC3<R CGAL_CTAG>::Direction_3
PlaneC3<R CGAL_CTAG>::orthogonal_direction() const
{
  return Direction_3(a(), b(), c());
}

template < class R >
typename PlaneC3<R CGAL_CTAG>::Vector_3
PlaneC3<R CGAL_CTAG>::base1() const
{
  if ( CGAL_NTS is_zero(a()) )  // parallel to x-axis
      return Vector_3(FT(1), FT(0), FT(0));

  if ( CGAL_NTS is_zero(b()) )  // parallel to y-axis
      return Vector_3(FT(0), FT(1), FT(0));

  if ( CGAL_NTS is_zero(c()) )  // parallel to z-axis
      return Vector_3(FT(0), FT(0), FT(1));

  return Vector_3(-b(), a(), FT(0)); // FIXME : construction
}

template < class R >
typename PlaneC3<R CGAL_CTAG>::Vector_3
PlaneC3<R CGAL_CTAG>::base2() const
{ // FIXME : construction
  if ( CGAL_NTS is_zero(a()) ) // parallel to x-axis already returned in base1
    {
      if ( CGAL_NTS is_zero(b()) )  // parallel to y-axis
          return Vector_3(FT(0), FT(1), FT(0));

      if ( CGAL_NTS is_zero(c()) ) // parallel to z-axis
          return Vector_3(FT(0), FT(0), FT(1));

      return Vector_3(FT(0), -b(), c());
    }
  if ( CGAL_NTS is_zero(b()) )
      return Vector_3(c(), FT(0), -a());

  if ( CGAL_NTS is_zero(c()) )
      return Vector_3(-b(), a(), FT(0));

  return Vector_3(FT(0), -c(), b());
}

template < class R >
typename PlaneC3<R CGAL_CTAG>::Point_3
PlaneC3<R CGAL_CTAG>::
to_plane_basis(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{ // FIXME : construction
  FT alpha, beta, gamma;

  solve(base1(), base2(), orthogonal_vector(), p - point(),
	alpha, beta, gamma);

  return Point_3(alpha, beta, gamma);
}

template < class R >
typename PlaneC3<R CGAL_CTAG>::Point_2
PlaneC3<R CGAL_CTAG>::
to_2d(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{ // FIXME : construction
  FT alpha, beta, gamma;

  solve(base1(), base2(), orthogonal_vector(), p - point(),
	alpha, beta, gamma);

  return Point_2(alpha, beta);
}

template < class R >
inline
typename PlaneC3<R CGAL_CTAG>::Point_3
PlaneC3<R CGAL_CTAG>::
to_3d(const typename PlaneC3<R CGAL_CTAG>::Point_2 &p) const
{ // FIXME : construction
  return point() + p.x() * base1() + p.y() * base2();
}

template < class R >
inline
typename PlaneC3<R CGAL_CTAG>::Line_3
PlaneC3<R CGAL_CTAG>::
perpendicular_line(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return Line_3(p, orthogonal_direction());
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>
PlaneC3<R CGAL_CTAG>::opposite() const
{ // FIXME : construction
  return PlaneC3<R CGAL_CTAG>(-a(), -b(), -c(), -d());
}

template < class R >
inline
Oriented_side
PlaneC3<R CGAL_CTAG>::
oriented_side(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return side_of_oriented_plane(*this, p);
}

#ifndef CGAL_NO_DEPRECATED_CODE
template < class R >
inline
bool
PlaneC3<R CGAL_CTAG>::
has_on_boundary(const  typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return has_on(p);
}

template < class R >
inline
bool
PlaneC3<R CGAL_CTAG>::
has_on_boundary(const  typename PlaneC3<R CGAL_CTAG>::Line_3 &l) const
{
  return has_on(l);
}
#endif

template < class R >
inline
bool
PlaneC3<R CGAL_CTAG>::
has_on(const  typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return oriented_side(p) == ON_ORIENTED_BOUNDARY;
}

template < class R >
inline
bool
PlaneC3<R CGAL_CTAG>::
has_on(const  typename PlaneC3<R CGAL_CTAG>::Line_3 &l) const
{
  return has_on(l.point())
     &&  has_on(l.point() + l.direction().to_vector());
}

template < class R >
inline
bool
PlaneC3<R CGAL_CTAG>::
has_on_positive_side(const  typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
inline
bool
PlaneC3<R CGAL_CTAG>::
has_on_negative_side(const  typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
PlaneC3<R CGAL_CTAG>::
is_degenerate() const
{ // FIXME : predicate
  return CGAL_NTS is_zero(a()) && CGAL_NTS is_zero(b()) &&
         CGAL_NTS is_zero(c());
}

#ifndef CGAL_NO_OSTREAM_INSERT_PLANEC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const PlaneC3<R CGAL_CTAG> &p)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << p.a() << ' ' << p.b() <<  ' ' << p.c() << ' ' << p.d();
    case IO::BINARY :
        write(os, p.a());
        write(os, p.b());
        write(os, p.c());
        write(os, p.d());
        return os;
        default:
            os << "PlaneC3(" << p.a() <<  ", " << p.b() <<   ", ";
            os << p.c() << ", " << p.d() <<")";
            return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_PLANEC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_PLANEC3
template < class R >
std::istream &
operator>>(std::istream &is, PlaneC3<R CGAL_CTAG> &p)
{
    typename R::FT a, b, c, d;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> a >> b >> c >> d;
        break;
    case IO::BINARY :
        read(is, a);
        read(is, b);
        read(is, c);
        read(is, d);
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
	p = PlaneC3<R CGAL_CTAG>(a, b, c, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_PLANEC3

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_PLANE_3_H
