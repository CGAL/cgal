// ======================================================================
//
// Copyright (c) 1998, 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.1-I-28 $
// release_date  : $CGAL_Date: 1999/10/12 $
//
// file          : include/CGAL/Cartesian/Plane_3.C
// package       : C3 (3.6.2)
// source        : include/CGAL/Cartesian/Plane_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

#include <CGAL/Cartesian/constructions_on_planes_3.h>
#include <CGAL/Cartesian/distance_computations_3.h>
#include <CGAL/Cartesian/predicates_on_planes_3.h>

#ifndef CGAL_CARTESIAN_PLANE_3_C
#define CGAL_CARTESIAN_PLANE_3_C

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Fourtuple<typename R::FT>*
PlaneC3<R CGAL_CTAG>::ptr() const
{
    return (_Fourtuple<FT>*)PTR;
}

template < class R >
inline
void
PlaneC3<R CGAL_CTAG>::
new_rep(const typename PlaneC3<R CGAL_CTAG>::FT &a,
        const typename PlaneC3<R CGAL_CTAG>::FT &b,
        const typename PlaneC3<R CGAL_CTAG>::FT &c,
        const typename PlaneC3<R CGAL_CTAG>::FT &d)
{
  PTR = new _Fourtuple<FT>(a, b, c, d);
}

template < class R >
inline
void
PlaneC3<R CGAL_CTAG>::new_rep(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p,
                              const typename PlaneC3<R CGAL_CTAG>::Point_3 &q,
                              const typename PlaneC3<R CGAL_CTAG>::Point_3 &r)
{
  PlaneC3<R CGAL_CTAG> h = plane_from_points(p,q,r);
  new_rep(h.a(), h.b(), h.c(), h.d());
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3()
{
  PTR = new _Fourtuple<FT>();
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const PlaneC3<R CGAL_CTAG> &p)
  : Handle(p)
{}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p,
        const typename PlaneC3<R CGAL_CTAG>::Point_3 &q,
        const typename PlaneC3<R CGAL_CTAG>::Point_3 &r)
{
  new_rep(p, q, r);
}

template < class R >
CGAL_KERNEL_INLINE
PlaneC3<R CGAL_CTAG>::
PlaneC3(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p,
        const typename PlaneC3<R CGAL_CTAG>::Direction_3 &d)
{
  PlaneC3<R CGAL_CTAG> h = plane_from_point_direction(p,d);
  new_rep(h.a(), h.b(), h.c(), h.d());
}

template < class R >
CGAL_KERNEL_INLINE
PlaneC3<R CGAL_CTAG>::
PlaneC3(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p,
        const typename PlaneC3<R CGAL_CTAG>::Vector_3 &v)
{
  FT a, b, c, d;
  plane_from_point_directionC3(p.x(),p.y(),p.z(),v.x(),v.y(),v.z(),a,b,c,d);
  new_rep(a, b, c, d);
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const typename PlaneC3<R CGAL_CTAG>::FT &a,
        const typename PlaneC3<R CGAL_CTAG>::FT &b,
        const typename PlaneC3<R CGAL_CTAG>::FT &c,
        const typename PlaneC3<R CGAL_CTAG>::FT &d)
{
  new_rep(a, b, c, d);
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const typename PlaneC3<R CGAL_CTAG>::Line_3 &l,
        const typename PlaneC3<R CGAL_CTAG>::Point_3 &p)
{
  new_rep(l.point(), l.point()+l.direction().to_vector(), p);
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const typename PlaneC3<R CGAL_CTAG>::Segment_3 &s,
        const typename PlaneC3<R CGAL_CTAG>::Point_3 &p)
{
  new_rep(s.start(), s.end(), p);
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::
PlaneC3(const typename PlaneC3<R CGAL_CTAG>::Ray_3 &r,
        const typename PlaneC3<R CGAL_CTAG>::Point_3 &p)
{
  new_rep(r.start(), r.second_point(), p);
}

template < class R >
inline
PlaneC3<R CGAL_CTAG>::~PlaneC3()
{}

template < class R >
inline
PlaneC3<R CGAL_CTAG> &PlaneC3<R CGAL_CTAG>::
operator=(const PlaneC3<R CGAL_CTAG> &p)
{
  Handle::operator=(p);
  return *this;
}

template < class R >
CGAL_KERNEL_INLINE
bool PlaneC3<R CGAL_CTAG>::
operator==(const PlaneC3<R CGAL_CTAG> &p) const
{
  return has_on_boundary(p.point()) &&
         (orthogonal_direction() == p.orthogonal_direction());

}

template < class R >
inline
bool PlaneC3<R CGAL_CTAG>::
operator!=(const PlaneC3<R CGAL_CTAG> &p) const
{
  return !(*this == p);
}

template < class R >
inline
long PlaneC3<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
inline
typename PlaneC3<R CGAL_CTAG>::FT
PlaneC3<R CGAL_CTAG>::a() const
{
  return ptr()->e0;
}

template < class R >
inline
typename PlaneC3<R CGAL_CTAG>::FT
PlaneC3<R CGAL_CTAG>::b() const
{
  return ptr()->e1;
}

template < class R >
inline
typename PlaneC3<R CGAL_CTAG>::FT
PlaneC3<R CGAL_CTAG>::c() const
{
  return ptr()->e2;
}

template < class R >
inline
typename PlaneC3<R CGAL_CTAG>::FT
PlaneC3<R CGAL_CTAG>::d() const
{
  return ptr()->e3;
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
  return CGAL::projection_plane(p, *this);
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
  if( a() == FT(0) )  // parallel to x-axis
      return Vector_3(FT(1), FT(0), FT(0));

  if( b() == FT(0) )  // parallel to y-axis
      return Vector_3(FT(0), FT(1), FT(0));

  if (c() == FT(0) )  // parallel to z-axis
      return Vector_3(FT(0), FT(0), FT(1));

  return Vector_3(-b(), a(), FT(0));
}

template < class R >
typename PlaneC3<R CGAL_CTAG>::Vector_3
PlaneC3<R CGAL_CTAG>::base2() const
{
  if ( a() == FT(0) ) // parallel to x-axis  x-axis already returned in base1
    {
      if (b() == FT(0) )  // parallel to y-axis
          return Vector_3(FT(0), FT(1), FT(0));

      if (c() == FT(0) ) // parallel to z-axis
          return Vector_3(FT(0), FT(0), FT(1));

      return Vector_3(FT(0), -b(), c());
    }
  if (b() == FT(0) )
      return Vector_3(c(), FT(0), -a());

  if (c() == FT(0) )
      return Vector_3(-b(), a(), FT(0));

  return Vector_3(FT(0), -c(), b());
}

template < class R >
typename PlaneC3<R CGAL_CTAG>::Point_3
PlaneC3<R CGAL_CTAG>::
to_plane_basis(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  const Vector_3 &e1 = base1();
  const Vector_3 &e2 = base2();
  FT alpha, beta, gamma;

  solve(e1, e2, orthogonal_vector(), p - point(), alpha, beta, gamma);

  return Point_3(alpha, beta, gamma);
}

template < class R >
typename PlaneC3<R CGAL_CTAG>::Point_2
PlaneC3<R CGAL_CTAG>::
to_2d(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  const Vector_3 &e1 = base1();
  const Vector_3 &e2 = base2();
  FT alpha, beta, gamma;

  solve(e1, e2, orthogonal_vector(), p - point(), alpha, beta, gamma);

  return Point_2(alpha, beta);
}

template < class R >
inline
typename PlaneC3<R CGAL_CTAG>::Point_3
PlaneC3<R CGAL_CTAG>::
to_3d(const typename PlaneC3<R CGAL_CTAG>::Point_2 &p) const
{
  Vector_3 e1 = base1(),
               e2 = base2();
  return point() + p.x() * e1 + p.y() * e2;
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
{
  return PlaneC3<R CGAL_CTAG>(-a(),-b(),-c(),-d());
}

template < class R >
PlaneC3<R CGAL_CTAG>
PlaneC3<R CGAL_CTAG>::
transform(const typename PlaneC3<R CGAL_CTAG>::Aff_transformation_3& t) const
{
  return PlaneC3<R CGAL_CTAG>( t.transform(point()), (t.is_even())
           ?   t.transpose().inverse().transform(orthogonal_direction())
           : - t.transpose().inverse().transform(orthogonal_direction()) );
}

template < class R >
inline
Oriented_side
PlaneC3<R CGAL_CTAG>::
oriented_side(const typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return side_of_oriented_plane(*this,p);
}

template < class R >
inline
bool
PlaneC3<R CGAL_CTAG>::
has_on_boundary(const  typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return oriented_side(p) == ON_ORIENTED_BOUNDARY;
}

template < class R >
inline
bool
PlaneC3<R CGAL_CTAG>::
has_on(const  typename PlaneC3<R CGAL_CTAG>::Point_3 &p) const
{
  return has_on_boundary(p);
}

template < class R >
inline
bool
PlaneC3<R CGAL_CTAG>::
has_on_boundary(const  typename PlaneC3<R CGAL_CTAG>::Line_3 &l) const
{
  return has_on_boundary(l.point())
     &&  has_on_boundary(l.point() + l.direction().to_vector());
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
{
  return (a() == FT(0)) && (b() == FT(0)) && (c() == FT(0));
}


#ifndef CGAL_NO_OSTREAM_INSERT_PLANEC3
template < class R >
std::ostream &operator<<(std::ostream &os, const PlaneC3<R CGAL_CTAG> &p)
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
std::istream &operator>>(std::istream &is, PlaneC3<R CGAL_CTAG> &p)
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
    p = PlaneC3<R CGAL_CTAG>(a, b, c, d);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_PLANEC3


CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_PLANE_3_C
