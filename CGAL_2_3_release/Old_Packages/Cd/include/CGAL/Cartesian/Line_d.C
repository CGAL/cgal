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
// file          : include/CGAL/Cartesian/Line_d.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Hervé Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_LINE_D_C
#define CGAL_CARTESIAN_LINE_D_C

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

#include <CGAL/Cartesian/Segment_d.h>
#include <CGAL/Cartesian/Ray_d.h>
#include <CGAL/Cartesian/Plane_d.h>
#include <CGAL/Cartesian/constructions_on_lines_d.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Twotuple< typename LineCd<R CGAL_CTAG>::Point_d > *
LineCd<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple<Point_d>*)PTR;
}

template < class R >
inline
void
LineCd<R CGAL_CTAG>::
new_rep(const typename LineCd<R CGAL_CTAG>::Point_d &p,
        const typename LineCd<R CGAL_CTAG>::Vector_d &v)
{
  CGAL_kernel_precondition(  p.dimension() == v.dimension() );
  // CGAL_kernel_precondition(  v != NULL_VECTOR );
  PTR = new _Twotuple< Point_d > (p, ORIGIN+v);
}

template < class R >
LineCd<R CGAL_CTAG>::LineCd()
{
  PTR = new _Twotuple< Point_d >();
}

template < class R >
LineCd<R CGAL_CTAG>::LineCd(const LineCd<R CGAL_CTAG>  &l)
  : Handle((Handle&)l)
{}

template < class R >
LineCd<R CGAL_CTAG>::
LineCd(const typename LineCd<R CGAL_CTAG>::Point_d &p,
       const typename LineCd<R CGAL_CTAG>::Point_d &q)
{
  new_rep(p, q-p);
}

template < class R >
LineCd<R CGAL_CTAG>::
LineCd(const typename LineCd<R CGAL_CTAG>::Segment_d &s)
{
  new_rep(s.start(), s.end() - s.start());
}

template < class R >
LineCd<R CGAL_CTAG>::
LineCd(const typename LineCd<R CGAL_CTAG>::Ray_d &r)
{
  new_rep(r.start(), r.point(1) - r.start());
}

template < class R >
inline
LineCd<R CGAL_CTAG>::
LineCd(const typename LineCd<R CGAL_CTAG>::Point_d &p,
       const typename LineCd<R CGAL_CTAG>::Direction_d &d)
{
  new_rep(p, d.to_vector());
}

template < class R >
LineCd<R CGAL_CTAG>::~LineCd()
{}


template < class R >
inline
LineCd<R CGAL_CTAG> &
LineCd<R CGAL_CTAG>::operator=(const LineCd<R CGAL_CTAG> &l)
{
  Handle::operator=(l);
  return *this;
}

template < class R >
inline
bool
LineCd<R CGAL_CTAG>::operator==(const LineCd<R CGAL_CTAG> &l) const
{
  if (ptr() == l.ptr()) return true; // identical
  return has_on(l.point()) && (direction() == l.direction());
}

template < class R >
inline
bool
LineCd<R CGAL_CTAG>::operator!=(const LineCd<R CGAL_CTAG> &l) const
{
  return !(*this == l);
}

template < class R >
inline
long
LineCd<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
inline
int
LineCd<R CGAL_CTAG>::dimension() const
{
  return ptr()->e0.dimension();
}

template < class R >
inline
typename LineCd<R CGAL_CTAG>::Point_d
LineCd<R CGAL_CTAG>::point() const
{
  return ptr()->e0;
}

template < class R >
inline
typename LineCd<R CGAL_CTAG>::Direction_d
LineCd<R CGAL_CTAG>::
direction() const
{
  return ((ptr()->e1) - ORIGIN).direction();
}


template < class R >
inline
typename LineCd<R CGAL_CTAG>::Point_d
LineCd<R CGAL_CTAG>::
point(int i) const
{
  return point_on_line(FT(i),*this);
}

template < class R >
inline
typename LineCd<R CGAL_CTAG>::Plane_d
LineCd<R CGAL_CTAG>::
perpendicular_plane(const typename LineCd<R CGAL_CTAG>::Point_d &p) const
{
  return Plane_d(p, direction().to_vector());
}

template < class R >
inline
LineCd<R CGAL_CTAG>
LineCd<R CGAL_CTAG>::
opposite() const
{
  return LineCd<R CGAL_CTAG>(point(), -direction());
}

template < class R >
inline
typename LineCd<R CGAL_CTAG>::Point_d
LineCd<R CGAL_CTAG>::
projection(const typename LineCd<R CGAL_CTAG>::Point_d &p) const
{
  return projection_line(p,*this);
}

template < class R >
inline
bool
LineCd<R CGAL_CTAG>::
has_on(const typename LineCd<R CGAL_CTAG>::Point_d &p) const
{
  return collinear(point(), point()+direction().to_vector(), p);
}


template < class R >
inline
bool
LineCd<R CGAL_CTAG>::is_degenerate() const
{
  return direction().is_degenerate();
}

template < class R >
inline
LineCd<R CGAL_CTAG>
LineCd<R CGAL_CTAG>::
transform(const typename LineCd<R CGAL_CTAG>::Aff_transformation_d &t) const
{
  return LineCd<R CGAL_CTAG>( t.transform(point()), t.transform(direction()));
}

#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_LINECD
template < class R >
std::ostream &operator<<(std::ostream &os, const LineCd<R CGAL_CTAG> &l)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << l.point(0) << ' ' << l.point(1);
    case IO::BINARY :
        return os << l.point(0) <<  l.point(1);
    default:
        return  os << "LineCd(" << l.point(0) << ", " << l.point(1) << ")";
    }
}
#endif // CGAL_CARTESIAN_NO_OSTREAM_INSERT_LINECD

#ifndef CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LINECD
template < class R >
std::istream &operator>>(std::istream &is, LineCd<R CGAL_CTAG> &l)
{
    typename LineCd<R CGAL_CTAG>::Point_d p, q;
    is >> p >> q;
    if (is)
        l = LineCd<R CGAL_CTAG>(p, q);
    return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LINECD

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_LINE_D_C
