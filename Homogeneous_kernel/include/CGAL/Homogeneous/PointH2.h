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

#ifndef CGAL_HOMOGENEOUS_POINT_2_H
#define CGAL_HOMOGENEOUS_POINT_2_H

#include <CGAL/Origin.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/array.h>
#include <CGAL/Kernel_d/Cartesian_const_iterator_d.h>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/utility.hpp>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointH2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::RT                   RT;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Direction_2          Direction_2;

  typedef boost::array<RT, 3>               Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  typedef Rational_traits<FT>  Rat_traits;

  Base base;

public:

  typedef FT Cartesian_coordinate_type;
  typedef const RT& Homogeneous_coordinate_type;
  typedef Cartesian_const_iterator_d<const RT*> Cartesian_const_iterator;
  typedef R_                                    R;

    PointH2() {}

    PointH2(const Origin &)
       : base(CGALi::make_array(RT(0), RT(0), RT(1))) {}

    template < typename Tx, typename Ty >
    PointH2(const Tx & x, const Ty & y,
            typename boost::enable_if< boost::mpl::and_<boost::is_convertible<Tx, RT>,
                                                        boost::is_convertible<Ty, RT> > >::type* = 0)
      : base(CGALi::make_array<RT>(x, y, RT(1))) {}

    PointH2(const FT& x, const FT& y)
      : base(CGALi::make_array<RT>(
             Rat_traits().numerator(x) * Rat_traits().denominator(y),
             Rat_traits().numerator(y) * Rat_traits().denominator(x),
             Rat_traits().denominator(x) * Rat_traits().denominator(y)))
    {
      CGAL_kernel_assertion(hw() > 0);
    }

    PointH2(const RT& hx, const RT& hy, const RT& hw)
      : base( hw >= RT(0) ? CGALi::make_array(hx, hy, hw)
                          : CGALi::make_array<RT>(-hx, -hy, -hw) ) {}

    bool    operator==( const PointH2<R>& p) const;
    bool    operator!=( const PointH2<R>& p) const;

    const RT & hx() const { return get(base)[0]; }
    const RT & hy() const { return get(base)[1]; }
    const RT & hw() const { return get(base)[2]; }

    FT      x()  const { return FT(hx()) / FT(hw()); }
    FT      y()  const { return FT(hy()) / FT(hw()); }

    FT      cartesian(int i)   const;
    FT      operator[](int i)  const;
    const RT & homogeneous(int i) const;

    Cartesian_const_iterator cartesian_begin() const
    {
      return make_cartesian_const_iterator_begin(get(base).begin(),
		                                 boost::prior(get(base).end()));
    }

    Cartesian_const_iterator cartesian_end() const
    {
      return make_cartesian_const_iterator_end(boost::prior(get(base).end()));
    }

    int     dimension() const;

    Direction_2 direction() const;
};

template < class R >
CGAL_KERNEL_INLINE
bool
PointH2<R>::operator==( const PointH2<R>& p) const
{ // FIXME : Predicate
  return (  (hx() * p.hw() == p.hx() * hw() )
          &&(hy() * p.hw() == p.hy() * hw() ) );
}

template < class R >
inline
bool
PointH2<R>::operator!=( const PointH2<R>& p) const
{ return !(*this == p); }

template < class R >
CGAL_KERNEL_INLINE
typename PointH2<R>::FT
PointH2<R>::cartesian(int i) const
{
  CGAL_kernel_precondition( (i==0 || i==1) );
  if (i==0)
      return x();
  return y();
}

template < class R >
CGAL_KERNEL_INLINE
const typename PointH2<R>::RT &
PointH2<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition( (i>=0) && (i<=2) );
  return get(base)[i];
}

template < class R >
inline
typename PointH2<R>::FT
PointH2<R>::operator[](int i) const
{ return cartesian(i); }


template < class R >
inline
int
PointH2<R>::dimension() const
{ return 2; }

template < class R >
CGAL_KERNEL_INLINE
typename PointH2<R>::Direction_2
PointH2<R>::direction() const
{ return typename PointH2<R>::Direction_2(*this); }

CGAL_END_NAMESPACE

#endif // CGAL_HOMOGENEOUS_POINT_2_H
