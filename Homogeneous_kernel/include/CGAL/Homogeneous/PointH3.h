// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_HOMOGENEOUS_POINT_3_H
#define CGAL_HOMOGENEOUS_POINT_3_H

#include <CGAL/Origin.h>
#include <CGAL/array.h>
#include <CGAL/Kernel/Cartesian_coordinate_iterator_3.h>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/logical.hpp>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointH3
{
   typedef typename R_::RT                   RT;
   typedef typename R_::FT                   FT;
   typedef typename R_::Vector_3             Vector_3;
   typedef typename R_::Point_3              Point_3;
   typedef typename R_::Direction_3          Direction_3;
   typedef typename R_::Aff_transformation_3 Aff_transformation_3;

   typedef boost::array<RT, 4>               Rep;
   typedef typename R_::template Handle<Rep>::type  Base;

   typedef Rational_traits<FT>  Rat_traits;

   Base base;

public:

   typedef Cartesian_coordinate_iterator_3<R_> Cartesian_const_iterator;
   typedef R_                 R;

  PointH3() {}

  PointH3(const Origin &)
    : base(CGALi::make_array(RT(0), RT(0), RT(0), RT(1))) {}

  template < typename Tx, typename Ty, typename Tz >
  PointH3(const Tx & x, const Ty & y, const Tz & z,
          typename boost::enable_if< boost::mpl::and_< boost::mpl::and_< boost::is_convertible<Tx, RT>,
                                                                         boost::is_convertible<Ty, RT> >,
                                                       boost::is_convertible<Tz, RT> > >::type* = 0)
    : base(CGALi::make_array<RT>(x, y, z, RT(1))) {}

  PointH3(const FT& x, const FT& y, const FT& z)
    : base(CGALi::make_array<RT>(
           Rat_traits().numerator(x) * Rat_traits().denominator(y)
                                     * Rat_traits().denominator(z),
           Rat_traits().numerator(y) * Rat_traits().denominator(x)
                                     * Rat_traits().denominator(z),
           Rat_traits().numerator(z) * Rat_traits().denominator(x)
                                     * Rat_traits().denominator(y),
           Rat_traits().denominator(x) * Rat_traits().denominator(y)
                                       * Rat_traits().denominator(z)))
  {
    CGAL_kernel_assertion(hw() > 0);
  }

  PointH3(const RT& x, const RT& y, const RT& z, const RT& w)
    : base( w < RT(0) ? CGALi::make_array<RT>(-x, -y, -z, -w)
                      : CGALi::make_array(x, y, z, w)) {}

  FT    x()  const;
  FT    y()  const;
  FT    z()  const;
  const RT & hx() const;
  const RT & hy() const;
  const RT & hz() const;
  const RT & hw() const;
  const RT & homogeneous(int i) const;
  FT    cartesian(int i) const;
  FT    operator[](int i) const;


  Cartesian_const_iterator cartesian_begin() const
  {
    return Cartesian_const_iterator(static_cast<const Point_3*>(this), 0);
  }

  Cartesian_const_iterator cartesian_end() const
  {
    return Cartesian_const_iterator(static_cast<const Point_3*>(this), 3);
  }

  int   dimension() const;

  Direction_3 direction() const;
  Point_3     transform( const Aff_transformation_3 & t) const;

  bool  operator==( const PointH3<R>& p) const;
  bool  operator!=( const PointH3<R>& p) const;
};


template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hx() const
{ return get(base)[0]; }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hy() const
{ return get(base)[1]; }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hz() const
{ return get(base)[2]; }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hw() const
{ return get(base)[3]; }

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::FT
PointH3<R>::x()  const
{ return FT(hx()) / FT(hw()); }

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::FT
PointH3<R>::y()  const
{ return FT(hy()) / FT(hw()); }

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::FT
PointH3<R>::z()  const
{ return FT(hz()) / FT(hw()); }

template < class R >
inline
int
PointH3<R>::dimension() const
{ return 3; }

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::FT
PointH3<R>::cartesian(int i) const
{
  CGAL_kernel_precondition(i == 0 || i == 1 || i == 2);
  return FT(get(base)[i]) / FT(hw());
}

template < class R >
CGAL_KERNEL_INLINE
const typename PointH3<R>::RT &
PointH3<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition(i == 0 || i == 1 || i == 2 || i == 3);
  return get(base)[i];
}

template < class R >
inline
typename PointH3<R>::FT
PointH3<R>::operator[](int i) const
{ return cartesian(i); }

template < class R >
inline
typename PointH3<R>::Direction_3
PointH3<R>::direction() const
{ return Direction_3(*this); }

template < class R >
CGAL_KERNEL_INLINE
bool
PointH3<R>::operator==( const PointH3<R> & p) const
{
  return ( (hx() * p.hw() == p.hx() * hw() )
         &&(hy() * p.hw() == p.hy() * hw() )
         &&(hz() * p.hw() == p.hz() * hw() ) );
}

template < class R >
inline
bool
PointH3<R>::operator!=( const PointH3<R> & p) const
{ return !(*this == p); }


template < class R >
inline
typename R::Point_3
PointH3<R>::transform(const typename PointH3<R>::Aff_transformation_3& t) const
{ return t.transform(static_cast<const Point_3&>(*this)); }

CGAL_END_NAMESPACE

#endif // CGAL_HOMOGENEOUS_POINT_3_H
