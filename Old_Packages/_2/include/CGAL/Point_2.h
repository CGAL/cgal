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

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif

#if defined CGAL_HOMOGENEOUS_H || defined CGAL_SIMPLE_HOMOGENEOUS_H
#include <CGAL/PointH2.h>
#endif

#if defined CGAL_CARTESIAN_H || defined CGAL_SIMPLE_CARTESIAN_H
#include <CGAL/Cartesian/Point_2.h>
#endif

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


  bool operator==(const CGAL::Point_2<R>& p) const
  {
    return RPoint_2::operator==(p);
  }

  bool operator!=(const CGAL::Point_2<R>& p) const
  {
    return !(*this == p);
  }

  RT hx() const
  {
    return RPoint_2::hx();
  }

  RT hy() const
  {
    return RPoint_2::hy();
  }

  RT hw() const
  {
    return RPoint_2::hw();
  }
  FT x() const
  {
    return RPoint_2::x();
  }

  FT y() const
  {
    return RPoint_2::y();
  }

  RT homogeneous(int i) const
  {
    return RPoint_2::homogeneous(i);
  }

  FT cartesian(int i) const
  {
    return RPoint_2::cartesian(i);
  }

  FT operator[](int i) const
  {
    return cartesian(i);
  }

  int dimension() const
  {
    return 2;
  }

  Bbox_2       bbox() const
  {
    return RPoint_2::bbox();
  }

  CGAL::Point_2<R> transform(const CGAL::Aff_transformation_2<R>& t) const
  {
    return RPoint_2::transform(t);
  }

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

template <class R>
inline
Point_2<R>
operator+(const Origin& o, const Vector_2<R>& v);

template <class R>
inline
Point_2<R>
operator-(const Origin& o, const Vector_2<R>& v);

template <class R>
inline
Vector_2<R>
operator-(const Point_2<R>& p, const Origin& );

CGAL_END_NAMESPACE

#include <CGAL/Vector_2.h>
#include <CGAL/point_vector_definitions_2.C>
#include <CGAL/Aff_transformation_2.h>

CGAL_BEGIN_NAMESPACE

template <class R>
inline
bool
operator==(const Origin& o, const Point_2<R>& p)
{ return p == o; }

template <class R>
inline
bool
operator!=(const Origin& o, const Point_2<R>& p)
{ return p != o; }

CGAL_END_NAMESPACE

#endif // CGAL_POINT_2_H
