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
// file          : Iso_rectangle_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_ISO_RECTANGLE_2_H
#define CGAL_ISO_RECTANGLE_2_H

CGAL_BEGIN_NAMESPACE

template <class R_>
class Iso_rectangle_2 : public R_::Kernel_base::Iso_rectangle_2
{
  typedef typename R_::RT                    RT;
  typedef typename R_::Point_2               Point_2;
  typedef typename R_::Kernel_base::Iso_rectangle_2  RIso_rectangle_2;
public:
  typedef  R_   R;

  Iso_rectangle_2()
    : RIso_rectangle_2() {}

  Iso_rectangle_2(const CGAL::Iso_rectangle_2<R> &r)
    : RIso_rectangle_2(r) {}

  Iso_rectangle_2(const RIso_rectangle_2& r)
    : RIso_rectangle_2(r) {}

  Iso_rectangle_2(const Point_2 &p, const Point_2 &q)
    : RIso_rectangle_2(p,q) {}

  Iso_rectangle_2(const Point_2 &left, const Point_2 &right,
                  const Point_2 &bottom, const Point_2 &top)
    : RIso_rectangle_2(left, right, bottom, top) {}

  Iso_rectangle_2(const RT& min_hx, const RT& min_hy, 
                  const RT& max_hx, const RT& max_hy)
    : RIso_rectangle_2(min_hx, min_hy, max_hx, max_hy) {}

  Iso_rectangle_2(const RT& min_hx, const RT& min_hy, 
                  const RT& max_hx, const RT& max_hy, const RT& hw)
    : RIso_rectangle_2(min_hx, min_hy, max_hx, max_hy, hw) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLE_2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Iso_rectangle_2<R> &r)
{
  typedef typename R::Kernel_base::Iso_rectangle_2  RIso_rectangle_2;
  return  os << (const RIso_rectangle_2&)r;
}
#endif // CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLE_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLE_2
template < class R >
std::istream &
operator>>(std::istream &is, Iso_rectangle_2<R> &r)
{
  typedef typename R::Kernel_base::Iso_rectangle_2  RIso_rectangle_2;
  is >> (RIso_rectangle_2&)r;
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLE_2

CGAL_END_NAMESPACE

#endif // CGAL_ISO_RECTANGLE_2_H
