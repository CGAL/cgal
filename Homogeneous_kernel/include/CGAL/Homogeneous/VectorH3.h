// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Stefan Schirra
 
#ifndef CGAL_HOMOGENEOUS_VECTOR_3_H
#define CGAL_HOMOGENEOUS_VECTOR_3_H

#include <CGAL/Origin.h>
#include <CGAL/array.h>
#include <CGAL/Kernel_d/Cartesian_const_iterator_d.h>

#include <boost/next_prior.hpp>

namespace CGAL {

template < class R_ >
class VectorH3
{
  typedef typename R_::RT                   RT;
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Segment_3            Segment_3;
  typedef typename R_::Ray_3                Ray_3;
  typedef typename R_::Line_3               Line_3;
  typedef typename R_::Direction_3          Direction_3;

  typedef cpp11::array<RT, 4>               Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  typedef Rational_traits<FT>               Rat_traits;

  Base base;

public:

  typedef Cartesian_const_iterator_d<typename Rep::const_iterator> Cartesian_const_iterator;

  typedef R_                 R;

  VectorH3() {}

  VectorH3(const Point_3& a, const Point_3& b)
  { *this = R().construct_vector_3_object()(a, b); }

  VectorH3(const Segment_3& s)
  { *this = R().construct_vector_3_object()(s); }

  VectorH3(const Ray_3& r)
  { *this = R().construct_vector_3_object()(r); }

  VectorH3(const Line_3& l)
  { *this = R().construct_vector_3_object()(l); }

  VectorH3(const Null_vector&)
    : base(CGAL::make_array(RT(0), RT(0), RT(0), RT(1))) {}

  template < typename Tx, typename Ty, typename Tz >
  VectorH3(const Tx & x, const Ty & y, const Tz & z,
           typename boost::enable_if< boost::mpl::and_< boost::mpl::and_< boost::is_convertible<Tx, RT>,
                                                                          boost::is_convertible<Ty, RT> >,
                                                        boost::is_convertible<Tz, RT> > >::type* = 0)
    : base(CGAL::make_array<RT>(x, y, z, RT(1))) {}

  VectorH3(const FT& x, const FT& y, const FT& z)
    : base(Rat_traits().denominator(x) * Rat_traits().denominator(y)
             * Rat_traits().denominator(z) >= 0 ?
               CGAL::make_array<RT>(
                 Rat_traits().numerator(x) * Rat_traits().denominator(y)
                                           * Rat_traits().denominator(z),
                 Rat_traits().numerator(y) * Rat_traits().denominator(x)
                                           * Rat_traits().denominator(z),
                 Rat_traits().numerator(z) * Rat_traits().denominator(x)
                                           * Rat_traits().denominator(y),
                 Rat_traits().denominator(x) * Rat_traits().denominator(y)
                                             * Rat_traits().denominator(z)) :
               CGAL::make_array<RT>(
                 - Rat_traits().numerator(x) * Rat_traits().denominator(y)
                                             * Rat_traits().denominator(z),
                 - Rat_traits().numerator(y) * Rat_traits().denominator(x)
                                             * Rat_traits().denominator(z),
                 - Rat_traits().numerator(z) * Rat_traits().denominator(x)
                                             * Rat_traits().denominator(y),
                 - Rat_traits().denominator(x) * Rat_traits().denominator(y)
                                               * Rat_traits().denominator(z)))
  {
    CGAL_kernel_assertion(hw() > 0);
  }

  VectorH3(const RT& x, const RT& y, const RT& z, const RT& w)
    : base( w >= RT(0) ? CGAL::make_array(x, y, z, w)
                       : CGAL::make_array<RT>(-x, -y, -z, -w) ) {}

  const RT & hx() const { return get_pointee_or_identity(base)[0]; }
  const RT & hy() const { return get_pointee_or_identity(base)[1]; }
  const RT & hz() const { return get_pointee_or_identity(base)[2]; }
  const RT & hw() const { return get_pointee_or_identity(base)[3]; }
  FT    x()  const { return FT(hx())/FT(hw()); }
  FT    y()  const { return FT(hy())/FT(hw()); }
  FT    z()  const { return FT(hz())/FT(hw()); }
  const RT & homogeneous(int i) const;
  FT    cartesian(int i) const;
  FT    operator[](int i) const;

  Cartesian_const_iterator cartesian_begin() const
  {
    return make_cartesian_const_iterator_begin(get_pointee_or_identity(base).begin(),
                                               boost::prior(get_pointee_or_identity(base).end()));
  }

  Cartesian_const_iterator cartesian_end() const
  {
    return make_cartesian_const_iterator_end(boost::prior(get_pointee_or_identity(base).end()));
  }

  int   dimension() const { return 3; };

  Direction_3 direction() const;

  Vector_3 operator-() const;

  bool  operator==( const VectorH3<R>& v) const;
  bool  operator!=( const VectorH3<R>& v) const;

  Vector_3 operator+( const VectorH3 &v) const;
  Vector_3 operator-( const VectorH3 &v) const;
  FT squared_length() const;
  Vector_3 operator/( const RT &f) const;
  Vector_3 operator/( const FT &f) const;
};


template < class R >
CGAL_KERNEL_INLINE
typename VectorH3<R>::FT
VectorH3<R>::cartesian(int i) const
{
  CGAL_kernel_precondition(i == 0 || i == 1 || i == 2);
  switch (i)
  {
      case 0:   return x();
      case 1:   return y();
  }
  return z();
}

template < class R >
CGAL_KERNEL_INLINE
const typename VectorH3<R>::RT &
VectorH3<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition(i == 0 || i == 1 || i == 2 || i == 3);
  return get_pointee_or_identity(base)[i];
}

template < class R >
inline
typename VectorH3<R>::Direction_3
VectorH3<R>::direction() const
{ return Direction_3(hx(), hy(), hz()); }

template < class R >
CGAL_KERNEL_INLINE
bool
VectorH3<R>::operator==( const VectorH3<R>& v) const
{
 return ( (hx() * v.hw() == v.hx() * hw() )
        &&(hy() * v.hw() == v.hy() * hw() )
        &&(hz() * v.hw() == v.hz() * hw() ) );
}

template < class R >
inline
bool
VectorH3<R>::operator!=( const VectorH3<R>& v) const
{ return !(*this == v); }

template < class R >
inline
typename VectorH3<R>::FT
VectorH3<R>::operator[](int i) const
{ return cartesian(i); }

template < class R >
CGAL_KERNEL_INLINE
typename VectorH3<R>::Vector_3
VectorH3<R>::operator-() const
{ return Vector_3( - hx(), - hy(), -hz(), hw() ); }

template <class R>
CGAL_KERNEL_INLINE
typename R::Vector_3
VectorH3<R>::operator+(const VectorH3<R>& v) const
{
  return typename R::Vector_3(hx()*v.hw() + v.hx()*hw(),
                              hy()*v.hw() + v.hy()*hw(),
                              hz()*v.hw() + v.hz()*hw(),
                              hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
typename R::Vector_3
VectorH3<R>::operator-(const VectorH3<R>& v) const
{
  return typename R::Vector_3(hx()*v.hw() - v.hx()*hw(),
                              hy()*v.hw() - v.hy()*hw(),
                              hz()*v.hw() - v.hz()*hw(),
                              hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
typename VectorH3<R>::FT
VectorH3<R>::squared_length() const
{
  typedef typename R::FT FT;
  return 
    FT( CGAL_NTS square(hx()) + 
	CGAL_NTS square(hy()) + 
	CGAL_NTS square(hz()) ) / 
    FT( CGAL_NTS square(hw()) );
}

template <class R>
CGAL_KERNEL_INLINE
typename R::Vector_3
VectorH3<R>::operator/(const typename VectorH3<R>::RT& f) const
{ return typename R::Vector_3( hx(), hy(), hz(), hw()*f ); }

template <class R>
CGAL_KERNEL_INLINE
typename R::Vector_3
VectorH3<R>::operator/(const typename VectorH3<R>::FT& f) const
{ return typename R::Vector_3(hx()*f.denominator(), hy()*f.denominator(),
		              hz()*f.denominator(), hw()*f.numerator() ); }

} //namespace CGAL

#endif // CGAL_HOMOGENEOUS_VECTOR_3_H
