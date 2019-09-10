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
 

#ifndef CGAL_HOMOGENEOUS_VECTOR_2_h
#define CGAL_HOMOGENEOUS_VECTOR_2_h

#include <CGAL/Origin.h>
#include <CGAL/array.h>
#include <CGAL/Kernel_d/Cartesian_const_iterator_d.h>
#include <CGAL/Handle_for.h>

#include <boost/next_prior.hpp>

namespace CGAL {

template < class R_ >
class VectorH2
{
  typedef VectorH2<R_>                      Self;
  typedef typename R_::FT                   FT;
  typedef typename R_::RT                   RT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Segment_2            Segment_2;
  typedef typename R_::Ray_2                Ray_2;
  typedef typename R_::Line_2               Line_2;
  typedef typename R_::Direction_2          Direction_2;
  typedef typename R_::Vector_2             Vector_2;

  typedef cpp11::array<RT, 3>               Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  typedef Rational_traits<FT>               Rat_traits;

  Base base;

public:

  typedef const FT Cartesian_coordinate_type;
  typedef const RT& Homogeneous_coordinate_type;
  typedef Cartesian_const_iterator_d<typename Rep::const_iterator> Cartesian_const_iterator;

  typedef R_                                    R;

   VectorH2() {}

   template < typename Tx, typename Ty >
   VectorH2(const Tx & x, const Ty & y,
            typename boost::enable_if< boost::mpl::and_<boost::is_convertible<Tx, RT>,
                                                        boost::is_convertible<Ty, RT> > >::type* = 0)
      : base(CGAL::make_array<RT>(x, y, RT(1))) {}

   VectorH2(const FT& x, const FT& y)
      : base(CGAL::make_array<RT>(
             Rat_traits().numerator(x) * Rat_traits().denominator(y),
             Rat_traits().numerator(y) * Rat_traits().denominator(x),
             Rat_traits().denominator(x) * Rat_traits().denominator(y)))
   {
     CGAL_kernel_assertion(hw() > 0);
   }

   VectorH2(const RT& x, const RT& y, const RT& w )
     : base( w >= RT(0) ? CGAL::make_array( x,  y,  w)
                        : CGAL::make_array<RT>(-x, -y, -w) ) {}

  const Self&
  rep() const
  {
    return static_cast<const Self& >(*this);
  }
  
   bool    operator==( const VectorH2<R>& v) const;
   bool    operator!=( const VectorH2<R>& v) const;
   bool    operator==( const Null_vector&) const;
   bool    operator!=( const Null_vector& v) const;

   const RT & hx() const { return CGAL::get_pointee_or_identity(base)[0]; };
   const RT & hy() const { return CGAL::get_pointee_or_identity(base)[1]; };
   const RT & hw() const { return CGAL::get_pointee_or_identity(base)[2]; };

   FT      x()  const { return FT(hx()) / FT(hw()); };
   FT      y()  const { return FT(hy()) / FT(hw()); };

   FT      cartesian(int i)   const;
   const RT & homogeneous(int i) const;
   FT      operator[](int i)  const;

   Cartesian_const_iterator cartesian_begin() const
   {
     return make_cartesian_const_iterator_begin(CGAL::get_pointee_or_identity(base).begin(),
                                                boost::prior(CGAL::get_pointee_or_identity(base).end()));
   }

   Cartesian_const_iterator cartesian_end() const
   {
     return make_cartesian_const_iterator_end(boost::prior(CGAL::get_pointee_or_identity(base).end()));
   }

   int     dimension() const;
   Direction_2 direction() const;
   Vector_2 perpendicular(const Orientation& o ) const;

  //   Vector_2 operator+(const VectorH2 &v) const;
   Vector_2 operator-(const VectorH2 &v) const;
   Vector_2 operator-() const;
   Vector_2 opposite() const;
   FT squared_length() const;
  //   Vector_2 operator/(const RT &f) const;
  //Vector_2 operator/(const FT &f) const;

// undocumented:
   VectorH2(const Direction_2 & dir)
      : base ( dir) {}

  VectorH2(const Point_2 & p)
     : base ( p) {}
};

template < class R >
inline
bool
VectorH2<R>::operator==( const Null_vector&) const
{ return (hx() == RT(0)) && (hy() == RT(0)); }

template < class R >
inline
bool
VectorH2<R>::operator!=( const Null_vector& v) const
{ return !(*this == v); }

template < class R >
CGAL_KERNEL_INLINE
bool
VectorH2<R>::operator==( const VectorH2<R>& v) const
{
  return (  (hx() * v.hw() == v.hx() * hw() )
          &&(hy() * v.hw() == v.hy() * hw() ) );
}

template < class R >
inline
bool
VectorH2<R>::operator!=( const VectorH2<R>& v) const
{ return !(*this == v); }  /* XXX */

template < class R >
CGAL_KERNEL_INLINE
typename VectorH2<R>::FT
VectorH2<R>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i==0 || i==1) );
  if (i==0)
      return x();
  return y();
}

template < class R >
CGAL_KERNEL_INLINE
const typename VectorH2<R>::RT &
VectorH2<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  return CGAL::get_pointee_or_identity(base)[i];
}

template < class R >
inline
typename VectorH2<R>::FT
VectorH2<R>::operator[](int i) const
{ return cartesian(i); }

template < class R >
inline
int
VectorH2<R>::dimension() const
{ return 2; }

template < class R >
CGAL_KERNEL_INLINE
typename VectorH2<R>::Direction_2
VectorH2<R>::direction() const
{ return Direction_2(hx(), hy()); }

template < class R >
inline
typename VectorH2<R>::Vector_2
VectorH2<R>::operator-() const
{ return VectorH2<R>(- hx(), - hy(), hw() ); }

template < class R >
inline
typename VectorH2<R>::Vector_2
VectorH2<R>::opposite() const
{ return VectorH2<R>(- hx(), - hy(), hw() ); }


template <class R>
CGAL_KERNEL_INLINE
typename VectorH2<R>::Vector_2
VectorH2<R>::operator-(const VectorH2<R>& v) const
{
  return VectorH2<R>( hx()*v.hw() - v.hx()*hw(),
                      hy()*v.hw() - v.hy()*hw(),
                      hw()*v.hw() );
}

template <class R>
CGAL_KERNEL_INLINE
typename VectorH2<R>::FT
VectorH2<R>::squared_length() const
{
  typedef typename R::FT FT;
  return 
    FT( CGAL_NTS square(hx()) + CGAL_NTS square(hy()) ) / 
    FT( CGAL_NTS square(hw()) );
}


template < class R >
CGAL_KERNEL_INLINE
typename R::Vector_2
VectorH2<R>::perpendicular(const Orientation& o) const
{
  CGAL_kernel_precondition(o != COLLINEAR);
  if (o == COUNTERCLOCKWISE)
      return typename R::Vector_2(-hy(), hx(), hw());
  else
      return typename R::Vector_2(hy(), -hx(), hw());
}

} //namespace CGAL

#endif // CGAL_HOMOGENEOUS_VECTOR_2_h
