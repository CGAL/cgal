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
// release_date  : 2000, July 30
// 
// source        : Iso_rectangle_2.fw
// file          : Iso_rectangle_2.h
// package       : _2 (3.6)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 3.6
// revision_date : 30 Jul 2000 
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_ISO_RECTANGLE_2_H
#define CGAL_ISO_RECTANGLE_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#ifndef CGAL_ISO_RECTANGLEH2_H
#include <CGAL/Iso_rectangleH2.h>
#endif // CGAL_ISO_RECTANGLEH2_H
#endif // CGAL_HOMOGENEOUS_H

#ifdef CGAL_CARTESIAN_H
#ifndef CGAL_ISO_RECTANGLEC2_H
#include <CGAL/Cartesian/Iso_rectangle_2.h>
#endif // CGAL_ISO_RECTANGLEC2_H
#endif // CGAL_CARTESIAN_H

#ifdef CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/SimpleCartesian/Iso_rectangleS2.h>
#endif // CGAL_SIMPLE_CARTESIAN_H


#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Iso_rectangle_2 : public R_::Iso_rectangle_2_base
{
public:
  typedef  R_   R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Iso_rectangle_2_base  RIso_rectangle_2;

  Iso_rectangle_2()
    : RIso_rectangle_2()
  {}

  Iso_rectangle_2(const CGAL::Iso_rectangle_2<R> &r)
    : RIso_rectangle_2(r)
  {}

  Iso_rectangle_2(const RIso_rectangle_2& r)
    : RIso_rectangle_2(r)
  {}

  Iso_rectangle_2(const CGAL::Point_2<R> &p,
                       const CGAL::Point_2<R> &q)
    : RIso_rectangle_2(p,q)
  {}


  bool
  operator==(const CGAL::Iso_rectangle_2<R> &r) const
  { return  vertex(0) == r.vertex(0) && vertex(2) == r.vertex(2); }

  bool
  operator!=(const CGAL::Iso_rectangle_2<R> &r) const
  { return !(*this == r); }


  CGAL::Point_2<R>
  min() const
  { return RIso_rectangle_2::min(); }

  CGAL::Point_2<R>
  max() const
  { return RIso_rectangle_2::max(); }

  FT
  xmin() const
  { return RIso_rectangle_2::xmin(); }

  FT
  ymin() const
  { return RIso_rectangle_2::ymin(); }

  FT
  xmax() const
  { return RIso_rectangle_2::xmax(); }

  FT
  ymax() const
  { return RIso_rectangle_2::ymax(); }

  CGAL::Point_2<R>
  vertex(int i) const
  { return RIso_rectangle_2::vertex(i); }

  CGAL::Point_2<R>
  operator[](int i) const
  { return vertex(i); }

  Bounded_side
  bounded_side(const CGAL::Point_2<R> &p) const
  { return RIso_rectangle_2::bounded_side(p); }

  bool
  has_on_boundary(const CGAL::Point_2<R> &p) const
  { return RIso_rectangle_2::has_on_boundary(p); }

  bool
  has_on_bounded_side(const CGAL::Point_2<R> &p) const
  { return RIso_rectangle_2::has_on_bounded_side(p); }

  bool
  has_on_unbounded_side(const CGAL::Point_2<R> &p) const
  { return RIso_rectangle_2::has_on_unbounded_side(p); }

  bool
  is_degenerate() const
  { return RIso_rectangle_2::is_degenerate(); }

  CGAL::Iso_rectangle_2<R>
  transform(const CGAL::Aff_transformation_2<R> &t) const
  { return  RIso_rectangle_2::transform(t); }

};

#ifndef NO_OSTREAM_INSERT_ISO_RECTANGLE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Iso_rectangle_2<R> &r)
{
  typedef typename R::Iso_rectangle_2_base  RIso_rectangle_2;
  return  os << (const RIso_rectangle_2&)r;
}
#endif // NO_OSTREAM_INSERT_ISO_RECTANGLE_2

#ifndef NO_ISTREAM_EXTRACT_ISO_RECTANGLE_2
template < class R >
std::istream &
operator>>(std::istream &is, Iso_rectangle_2<R> &r)
{
  typedef typename R::Iso_rectangle_2_base  RIso_rectangle_2;
  is >> (RIso_rectangle_2&)r;
  return is;
}
#endif // NO_ISTREAM_EXTRACT_ISO_RECTANGLE_2

CGAL_END_NAMESPACE


#endif // CGAL_ISO_RECTANGLE_2_H
