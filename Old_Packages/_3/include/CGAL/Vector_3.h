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
// file          : Vector_3.h
// package       : _3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_VECTOR_3_H
#define CGAL_VECTOR_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#include <CGAL/Point_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Aff_transformation_3.h>

#ifdef VECTOR_WRAPPER
#include <CGAL/Vector_3_rft_wrapper.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class T> class Quotient;
template <class R_>
class Vector_3 : public R_::Vector_3_base
{
public:
  typedef          R_                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Vector_3_base  RVector_3;

friend CGAL_FRIEND_INLINE
       CGAL::Vector_3<R>
       point_to_vector_conversion CGAL_NULL_TMPL_ARGS
                                       (const CGAL::Point_3<R>& p);
/*
friend CGAL::Vector_3<R>
       CGAL::Direction_3<R>::vector() const;
*/

  Vector_3()
  {}
  Vector_3(const CGAL::Vector_3<R>& v)
    : RVector_3( (const RVector_3& )v )
  {}
  Vector_3(const RVector_3&  v) : RVector_3(v)
  {}
  Vector_3(const Null_vector& v) : RVector_3(v)
  {}
  Vector_3(const RT& x, const RT& y, const RT& z)
    : RVector_3(x, y, z)
  {}
  Vector_3(const RT& x, const RT& y, const RT& z, const RT& w)
    : RVector_3(x, y, z, w)
  {}

  bool operator==(const CGAL::Vector_3<R>& v) const
  { return RVector_3::operator==(v); }
  bool operator!=(const CGAL::Vector_3<R>& v) const
  { return !(*this == v); }
  bool operator==(const Null_vector& v) const
  { return RVector_3::operator==(v); }
  bool operator!=(const Null_vector& v) const
  { return !(*this == v); }
  RT hx() const
  { return RVector_3::hx(); }
  RT hy() const
  { return RVector_3::hy(); }
  RT hz() const
  { return RVector_3::hz(); }
  RT hw() const
  { return RVector_3::hw(); }
  FT x() const
  { return RVector_3::x(); }
  FT y() const
  { return RVector_3::y(); }
  FT z() const
  { return RVector_3::z(); }
  RT homogeneous(int i) const
  { return RVector_3::homogeneous(i); }
  FT cartesian(int i) const
  { return RVector_3::cartesian(i); }
  FT operator[](int i) const
  { return cartesian(i); }
  int dimension() const
  { return 3; }
  CGAL::Vector_3<R> operator+(const CGAL::Vector_3<R>& w) const
  { return (const RVector_3& )(*this) + (const RVector_3& )(w); }
  CGAL::Vector_3<R> operator-(const CGAL::Vector_3<R>& w) const
  { return (const RVector_3& )(*this) - (const RVector_3& )(w); }
  CGAL::Vector_3<R> operator-() const
  { return RVector_3::operator-(); }
  FT operator*(const CGAL::Vector_3<R>& w) const
  { return (const RVector_3& )(*this) * (const RVector_3& )(w); }

#ifndef VECTOR_WRAPPER
  CGAL::Vector_3<R> operator*(const RT& c) const
  { return c * (const RVector_3& )(*this) ; }
  CGAL::Vector_3<R> operator*(const Quotient<RT>& q) const
  {
    return (q.numerator() * (const RVector_3& )(*this)) /
            q.denominator();
  }
  CGAL::Vector_3<R> operator/(const Quotient<RT>& q) const
  {
    return (q.denominator() * (const RVector_3& )(*this)) /
            q.numerator();
  }
#endif // VECTOR_WRAPPER

  CGAL::Vector_3<R> operator/(const RT& c) const
  { return (const RVector_3& )(*this) / c; }
  CGAL::Direction_3<R> direction() const
  { return RVector_3::direction(); }
  CGAL::Vector_3<R> transform(const CGAL::Aff_transformation_3<R>& t) const
  { return RVector_3::transform(t); }

private:
  Vector_3(const CGAL::Point_3<R>& p) : RVector_3(p)
  {}
  Vector_3(const CGAL::Direction_3<R>& d) : RVector_3(d)
  {}
};

template < class R >
No_number_tag number_type_tag(const Vector_3<R>& )
{
  return No_number_tag();
}

#ifndef CGAL_NO_OSTREAM_INSERT_VECTOR_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Vector_3<R>& v)
{
  typedef typename  R::Vector_3_base  RVector_3;
  return os << (const RVector_3& )v;
}
#endif // CGAL_NO_OSTREAM_INSERT_VECTOR_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_VECTOR_3
template < class R >
std::istream&
operator>>(std::istream& is, Vector_3<R>& p)
{
  typedef typename  R::Vector_3_base  RVector_3;
  return is >> (RVector_3& )p;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_VECTOR_3

template<class R>
inline
Vector_3<R>
cross_product(const Vector_3<R>& v, const Vector_3<R>& w)
{
  typedef typename  R::Vector_3_base  RVector_3;
  return cross_product((const RVector_3& )v,(const RVector_3& )w);
}

CGAL_END_NAMESPACE

#endif // CGAL_VECTOR_3_H
