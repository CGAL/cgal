// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : Vector_3.fw
// file          : Vector_3.h
// revision      : 2.4
// revision_date : 24 Aug 1999 
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_VECTOR_3_H
#define CGAL_VECTOR_3_H

#ifndef CGAL_REP_CLASS_DEFINED
#error  no representation class defined
#endif  // CGAL_REP_CLASS_DEFINED

#ifndef CGAL_POINT_3_H
#include <CGAL/Point_3.h>
#endif // CGAL_POINT_3_H


#ifndef CGAL_DIRECTION_3_H
#include <CGAL/Direction_3.h>
#endif // CGAL_DIRECTION_3_H

#ifdef VECTOR_WRAPPER
#ifndef VECTOR_3_RFT_WRAPPER_H
#include <CGAL/Vector_3_rft_wrapper.h>
#endif // VECTOR_3_RFT_WRAPPER_H
#endif // VECTOR_WRAPPER

CGAL_BEGIN_NAMESPACE

template <class T> class Quotient;
template <class _R>
class Vector_3 : public _R::Vector_3_base
{
public:
  typedef          _R                       R;
  typedef typename R::RT                    RT;
  typedef typename R::FT                    FT;
  typedef typename R::Vector_3_base  RVector_3;

friend CGAL_FRIEND_INLINE
       Vector_3<R>
       point_to_vector_conversion CGAL_NULL_TMPL_ARGS 
                                       (const Point_3<R>& p);

friend Vector_3<R>
       Direction_3<R>::vector() const;

  Vector_3()
  {}
  Vector_3(const Vector_3<R>& v)
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

  Vector_3<R>& operator=(const Vector_3<R>& v)
  {
    RVector_3::operator=(v);
    return *this;
  }
  bool operator==(const Vector_3<R>& v) const
  { return RVector_3::operator==(v); }
  bool operator!=(const Vector_3<R>& v) const
  { return !(*this == v); }
  bool operator==(const Null_vector& v) const
  { return RVector_3::operator==(v); }
  bool operator!=(const Null_vector& v) const
  { return !(*this == v); }
  int id() const    /* XXX */
  { return (int) PTR ; }
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
  Vector_3<R> operator+(const Vector_3<R>& w) const
  { return (const RVector_3& )(*this) + (const RVector_3& )(w); }
  Vector_3<R> operator-(const Vector_3<R>& w) const
  { return (const RVector_3& )(*this) - (const RVector_3& )(w); }
  Vector_3<R> operator-() const
  { return RVector_3::operator-(); }
  FT operator*(const Vector_3<R>& w) const
  { return (const RVector_3& )(*this) * (const RVector_3& )(w); }

#ifndef VECTOR_WRAPPER
  Vector_3<R> operator*(const RT& c) const
  { return c * (const RVector_3& )(*this) ; }
  Vector_3<R> operator*(const Quotient<RT>& q) const
  {
    return (q.numerator() * (const RVector_3& )(*this)) /
            q.denominator();
  }
  Vector_3<R> operator/(const Quotient<RT>& q) const
  {
    return (q.denominator() * (const RVector_3& )(*this)) /
            q.numerator();
  }
#endif // VECTOR_WRAPPER

  Vector_3<R> operator/(const RT& c) const
  { return (const RVector_3& )(*this) / c; }
  Direction_3<R> direction() const
  { return RVector_3::direction(); }
  Vector_3<R> transform(const Aff_transformation_3<R>& t) const
  { return RVector_3::transform(t); }

private:
  Vector_3(const Point_3<R>& p) : RVector_3(p)
  {}
  Vector_3(const Direction_3<R>& d) : RVector_3(d)
  {}
};

template < class R >
No_number_tag number_type_tag(const Vector_3<R>& )
{
  return No_number_tag();
}

#ifdef VECTOR_WRAPPER
template <class T, class R>
_Vector_3_rft_wrapper<R>
multiply(const Quotient<T>& q,
              const Vector_3<R>& w,
              const Quotient_tag&)
{
  typedef typename  R::Vector_3_base  RVector_3;
  return _Vector_3_rft_wrapper<R>(
                 Vector_3<R>((q.numerator() * (const RVector_3& )(w))
                                  / q.denominator()));
}

template < class R >
_Vector_3_rft_wrapper<R>
multiply(const Vector_3<R>& v,
              const Vector_3<R>& w,
              const No_number_tag&)
{
  typedef typename  R::Vector_3_base  RVector_3;
  return _Vector_3_rft_wrapper<R>((const RVector_3& )(v)
                                     * (const RVector_3& )(w));
}

template < class T, class R >
_Vector_3_rft_wrapper<R>
multiply(const T& n,
              const Vector_3<R>& w,
              const Number_tag&)
{
  typedef typename  R::Vector_3_base  RVector_3;
  typedef typename  R::RT             RT;
  return _Vector_3_rft_wrapper<R>(
                 Vector_3<R>(RT(n) * (const RVector_3& )(w)));
}

template <class T, class R>
_Vector_3_rft_wrapper<R>
operator*(const T& t, const Vector_3<R>& w)
{
  return multiply(t, w, number_type_tag(t));
}
#endif // VECTOR_WRAPPER

#ifndef NO_OSTREAM_INSERT_VECTOR_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const Vector_3<R>& v)
{
  typedef typename  R::Vector_3_base  RVector_3;
  return os << (const RVector_3& )v;
}
#endif // NO_OSTREAM_INSERT_VECTOR_3

#ifndef NO_ISTREAM_EXTRACT_VECTOR_3
template < class R >
std::istream&
operator>>(std::istream& is, Vector_3<R>& p)
{
  typedef typename  R::Vector_3_base  RVector_3;
  return is >> (RVector_3& )p;
}
#endif // NO_ISTREAM_EXTRACT_VECTOR_3


template<class R>
inline
Vector_3<R>
cross_product(const Vector_3<R>& v,const Vector_3<R>& w)
{
  typedef typename  R::Vector_3_base  RVector_3;
  return cross_product((const RVector_3& )v,(const RVector_3& )w);
}

CGAL_END_NAMESPACE


#endif // CGAL_VECTOR_3_H
