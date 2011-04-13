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
// release_date  : 
// 
// file          : Iso_cuboid_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_ISO_CUBOID_3_H
#define CGAL_ISO_CUBOID_3_H

#include <CGAL/Point_3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Iso_cuboid_3 : public R_::Iso_cuboid_3_base
{
public:
  typedef          R_                    R;
  typedef typename R::RT                 RT;
  typedef typename R::FT                 FT;
  typedef typename R::Iso_cuboid_3_base  RIso_cuboid_3;

  Iso_cuboid_3() : RIso_cuboid_3()
  {}

  Iso_cuboid_3(const CGAL::Iso_cuboid_3<R>& r) : RIso_cuboid_3(r)
  {}

  Iso_cuboid_3(const RIso_cuboid_3&  r) : RIso_cuboid_3(r)
  {}

  Iso_cuboid_3(const CGAL::Point_3<R>& p, const CGAL::Point_3<R>& q)
   : RIso_cuboid_3(p,q)
  {}

  Iso_cuboid_3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz, 
               const RT& hw)
   : RIso_cuboid_3(min_hx, min_hy, min_hz, max_hx, max_hy, max_hz, hw)
  {}

  Iso_cuboid_3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz)
   : RIso_cuboid_3(min_hx, min_hy, min_hz, max_hx, max_hy, max_hz)
  {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_ISO_CUBOID_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Iso_cuboid_3<R>& r)
{
  typedef typename  R::Iso_cuboid_3_base  RIso_cuboid_3;
  return  os << (const RIso_cuboid_3& )r; }
#endif // CGAL_NO_OSTREAM_INSERT_ISO_CUBOID_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_ISO_CUBOID_3
template < class R >
std::istream&
operator>>(std::istream& is, Iso_cuboid_3<R>& r)
{
  typedef typename  R::Iso_cuboid_3_base  RIso_cuboid_3;
  is >> (RIso_cuboid_3& )r;
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISO_CUBOID_3

CGAL_END_NAMESPACE

#endif // CGAL_ISO_CUBOID_3_H
