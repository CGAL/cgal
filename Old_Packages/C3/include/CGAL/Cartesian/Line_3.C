// ==========================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// --------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Line_3.C
// source        : include/CGAL/Cartesian/Line_3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ==========================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

#ifndef CGAL_CARTESIAN_LINE_3_C
#define CGAL_CARTESIAN_LINE_3_C

#ifndef CGAL_CARTESIAN_SEGMENT_3_H
#include <CGAL/Cartesian/Segment_3.h>
#endif // CGAL_CARTESIAN_SEGMENT_3_H
#ifndef CGAL_CARTESIAN_RAY_3_H
#include <CGAL/Cartesian/Ray_3.h>
#endif // CGAL_CARTESIAN_RAY_3_H
#ifndef CGAL_CARTESIAN_PLANE_3_H
#include <CGAL/Cartesian/Plane_3.h>
#endif // CGAL_CARTESIAN_PLANE_3_H
#ifndef CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_3_H
#include <CGAL/Cartesian/constructions_on_lines_3.h>
#endif // CGAL_CARTESIAN_CONSTRUCTIONS_ON_LINES_3_H

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Twotuple< typename LineC3<R CGAL_CTAG>::Point_3 > *
LineC3<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple<Point_3>*)PTR;
}

template < class R >
inline
void
LineC3<R CGAL_CTAG>::
new_rep(const typename LineC3<R CGAL_CTAG>::Point_3 &p,
        const typename LineC3<R CGAL_CTAG>::Vector_3 &v)
{
  // CGAL_kernel_precondition(  v != NULL_VECTOR );
  PTR = new _Twotuple< Point_3 > (p, ORIGIN+v);
}

template < class R >
LineC3<R CGAL_CTAG>::LineC3()
{
  PTR = new _Twotuple<FT>();
}

template < class R >
LineC3<R CGAL_CTAG>::LineC3(const LineC3<R CGAL_CTAG>  &l)
  : Handle((Handle&)l)
{}

template < class R >
LineC3<R CGAL_CTAG>::
LineC3(const typename LineC3<R CGAL_CTAG>::Point_3 &p,
       const typename LineC3<R CGAL_CTAG>::Point_3 &q)
{
  new_rep(p, q-p);
}

template < class R >
LineC3<R CGAL_CTAG>::
LineC3(const typename LineC3<R CGAL_CTAG>::Segment_3 &s)
{
  new_rep(s.start(), s.end() - s.start());
}

template < class R >
LineC3<R CGAL_CTAG>::
LineC3(const typename LineC3<R CGAL_CTAG>::Ray_3 &r)
{
  new_rep(r.start(), r.point(1) - r.start());
}

template < class R >
LineC3<R CGAL_CTAG>::
LineC3(const typename LineC3<R CGAL_CTAG>::Point_3 &p,
       const typename LineC3<R CGAL_CTAG>::Direction_3 &d)
{
  new_rep(p, d.to_vector());
}

template < class R >
LineC3<R CGAL_CTAG>::~LineC3()
{}


template < class R >
LineC3<R CGAL_CTAG> &
LineC3<R CGAL_CTAG>::operator=(const LineC3<R CGAL_CTAG> &l)
{
  Handle::operator=(l);
  return *this;
}

template < class R >
bool
LineC3<R CGAL_CTAG>::operator==(const LineC3<R CGAL_CTAG> &l) const
{
  return has_on(l.point()) && (direction() == l.direction());
}

template < class R >
inline
bool
LineC3<R CGAL_CTAG>::operator!=(const LineC3<R CGAL_CTAG> &l) const
{
  return !(*this == l);
}

template < class R >
long
LineC3<R CGAL_CTAG>::id() const
{
  return (long) PTR ;
}

template < class R >
typename LineC3<R CGAL_CTAG>::Point_3
LineC3<R CGAL_CTAG>::point() const
{
  return ptr()->e0;
}

template < class R >
typename LineC3<R CGAL_CTAG>::Direction_3
LineC3<R CGAL_CTAG>::direction() const
{
  return ((ptr()->e1) - ORIGIN).direction();
}


template < class R >
typename LineC3<R CGAL_CTAG>::Point_3
LineC3<R CGAL_CTAG>::point(int i) const
{
  return point_on_line(i,*this);
}

template < class R >
typename LineC3<R CGAL_CTAG>::Plane_3
LineC3<R CGAL_CTAG>::
perpendicular_plane(const typename LineC3<R CGAL_CTAG>::Point_3 &p) const
{
  return Plane_3(p, direction().to_vector());
}

template < class R >
LineC3<R CGAL_CTAG>
LineC3<R CGAL_CTAG>::opposite() const
{
  return LineC3<R CGAL_CTAG>(point(), -direction());
}

template < class R >
typename LineC3<R CGAL_CTAG>::Point_3
LineC3<R CGAL_CTAG>::projection(const typename LineC3<R CGAL_CTAG>::Point_3 &p) const
{
  return projection_line(p,*this);
}

template < class R >
bool
LineC3<R CGAL_CTAG>::
has_on(const typename LineC3<R CGAL_CTAG>::Point_3 &p) const
{
  return collinear(point(), point()+direction().to_vector(), p);
}


template < class R >
bool
LineC3<R CGAL_CTAG>::is_degenerate() const
{
  return direction() == Direction_3(0,0,0);
}

template < class R >
inline
LineC3<R CGAL_CTAG>
LineC3<R CGAL_CTAG>::
transform(const typename LineC3<R CGAL_CTAG>::Aff_transformation_3 &t) const
{
  return LineC3<R CGAL_CTAG>( t.transform(point()), t.transform(direction()));
}


#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_LINEC3
template < class R >
std::ostream &operator<<(std::ostream &os, const LineC3<R CGAL_CTAG> &l)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << l.point(0) << ' ' << l.point(1);
    case IO::BINARY :
        return os << l.point(0) <<  l.point(1);
    default:
        return  os << "LineC3(" << l.point(0) << ", " << l.point(1) << ")";
    }
}
#endif // CGAL_CARTESIAN_NO_OSTREAM_INSERT_LINEC3

#ifndef CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LINEC3
template < class R >
std::istream &operator>>(std::istream &is, LineC3<R CGAL_CTAG> &l)
{
    typename LineC3<R CGAL_CTAG>::Point_3 p, q;
    is >> p >> q;
    l = LineC3<R CGAL_CTAG>(p, q);
    return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LINEC3

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_LINE_3_C
