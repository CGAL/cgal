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
// file          : include/CGAL/Cartesian/Line_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_LINE_3_H
#define CGAL_CARTESIAN_LINE_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Cartesian/Line_rep_3.h>
#include <CGAL/Cartesian/point_constructions_3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class LineC3 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Line_handle_3
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Line_handle_3             Line_handle_3_;
  typedef typename Line_handle_3_::element_type  Line_ref_3;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef LineC3<R CGAL_CTAG>                   Self;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Vector_3                  Vector_3;
  typedef typename R::Direction_3               Direction_3;
  typedef typename R::Plane_3                   Plane_3;
  typedef typename R::Ray_3                     Ray_3;
  typedef typename R::Segment_3                 Segment_3;
  typedef typename R::Aff_transformation_3      Aff_transformation_3;
#else
  typedef LineC3<R>                             Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Vector_3_base             Vector_3;
  typedef typename R::Direction_3_base          Direction_3;
  typedef typename R::Plane_3_base              Plane_3;
  typedef typename R::Ray_3_base                Ray_3;
  typedef typename R::Segment_3_base            Segment_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  LineC3()
    : Line_handle_3_(Line_ref_3()) {}

  LineC3(const Point_3 &p, const Point_3 &q) // FIXME : construction
    : Line_handle_3_(Line_ref_3(p, (q-p).direction())) {}

  LineC3(const Segment_3 &s) // FIXME : construction
    : Line_handle_3_(Line_ref_3(s.start(),
		               (s.end() - s.start()).direction())) {}

  LineC3(const Ray_3 &r) // FIXME : construction
    : Line_handle_3_(Line_ref_3(r.start(),
	                       (r.point(1) - r.start()).direction())) {}

  LineC3(const Point_3 &p, const Direction_3 &d)
    : Line_handle_3_(Line_ref_3(p, d)) {}

  bool        operator==(const Self &l) const;
  bool        operator!=(const Self &l) const;

  Plane_3     perpendicular_plane(const Point_3 &p) const;
  Self        opposite() const;

  Point_3     point() const
  {
      return Ptr()->basepoint;
  }
  Direction_3 direction() const
  {
      return Ptr()->direction;
  }

  Point_3     point(int i) const;

  Point_3     projection(const Point_3 &p) const;

  bool        has_on(const Point_3 &p) const;
  bool        is_degenerate() const;

  Self        transform(const Aff_transformation_3 &t) const
  {
    return Self(t.transform(point()), t.transform(direction()));
  }
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
inline
bool
LineC3<R CGAL_CTAG>::operator==(const LineC3<R CGAL_CTAG> &l) const
{
  if (identical(l))
      return true;
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
inline
typename LineC3<R CGAL_CTAG>::Point_3
LineC3<R CGAL_CTAG>::point(int i) const
{
  return point_on_line(i, *this);
}

template < class R >
inline
typename LineC3<R CGAL_CTAG>::Plane_3
LineC3<R CGAL_CTAG>::
perpendicular_plane(const typename LineC3<R CGAL_CTAG>::Point_3 &p) const
{
  return Plane_3(p, direction().to_vector());
}

template < class R >
inline
LineC3<R CGAL_CTAG>
LineC3<R CGAL_CTAG>::opposite() const
{
  return LineC3<R CGAL_CTAG>(point(), -direction());
}

template < class R >
inline
typename LineC3<R CGAL_CTAG>::Point_3
LineC3<R CGAL_CTAG>::
projection(const typename LineC3<R CGAL_CTAG>::Point_3 &p) const
{
  return projection_line(p, *this);
}

template < class R >
inline
bool
LineC3<R CGAL_CTAG>::
has_on(const typename LineC3<R CGAL_CTAG>::Point_3 &p) const
{
  return collinear(point(), point()+direction().to_vector(), p);
}

template < class R >
inline
bool
LineC3<R CGAL_CTAG>::is_degenerate() const
{ // FIXME : predicate
  return direction() == Direction_3(0,0,0);
}

#ifndef CGAL_CARTESIAN_NO_OSTREAM_INSERT_LINEC3
template < class R >
std::ostream &
operator<<(std::ostream &os, const LineC3<R CGAL_CTAG> &l)
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
std::istream &
operator>>(std::istream &is, LineC3<R CGAL_CTAG> &l)
{
    typename LineC3<R CGAL_CTAG>::Point_3 p, q;
    is >> p >> q;
    if (is)
	l = LineC3<R CGAL_CTAG>(p, q);
    return is;
}
#endif // CGAL_CARTESIAN_NO_ISTREAM_EXTRACT_LINEC3

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_LINE_3_H
