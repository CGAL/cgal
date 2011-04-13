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
// file          : Point_2.h
// package       : _2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_POINT_2_H
#define CGAL_POINT_2_H

#include <CGAL/point_vector_declarations_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Point_2 : public R_::Point_2_base
{
public:
  typedef  R_   R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Point_2_base  RPoint_2;
  typedef typename R::Vector_2_base  RVector_2;


friend  CGAL_FRIEND_INLINE
        CGAL::Point_2<R>
        CGAL_SCOPE vector_to_point_conversion CGAL_NULL_TMPL_ARGS
                                         (const CGAL::Vector_2<R>& v);

  Point_2()
  {}

  Point_2(const Origin& o)
    : RPoint_2(o)
  {}

  Point_2(const CGAL::Point_2<R>& p)
    : RPoint_2(static_cast<const RPoint_2&>(p))
  {}

  Point_2(const RPoint_2& p)
    : RPoint_2(p)
  {}

  Point_2(const RT& hx, const RT& hy)
    : RPoint_2(hx, hy)
  {}

  Point_2(const RT& hx, const RT& hy, const RT& hw)
    : RPoint_2(hx, hy, hw)
  {}

private:

  Point_2(const RVector_2& v)
    : RPoint_2(v)
  {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_POINT_2
template < class R >
std::ostream&
operator<<(std::ostream& os, const Point_2<R>& p)
{
  typedef typename  R::Point_2_base  RPoint_2;
  return os << static_cast<const RPoint_2&>(p);
}
#endif // CGAL_NO_OSTREAM_INSERT_POINT_2

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINT_2
template < class R >
std::istream&
operator>>(std::istream& is, Point_2<R>& p)
{
  typedef typename  R::Point_2_base  RPoint_2;
  return is >> static_cast<RPoint_2&>(p);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINT_2

CGAL_END_NAMESPACE

#include <CGAL/Vector_2.h>
#include <CGAL/Aff_transformation_2.h>

#endif // CGAL_POINT_2_H
