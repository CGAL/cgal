// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_HOMOGENEOUS_POINT_2_H
#define CGAL_HOMOGENEOUS_POINT_2_H

#include <CGAL/Origin.h>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/logical.hpp>

namespace CGAL {

template < class R_ >
class PointH2
{
  typedef typename R_::FT                   FT;
  typedef typename R_::RT                   RT;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Direction_2          Direction_2;

  typedef Rational_traits<FT>  Rat_traits;

  // Reference-counting is handled in Vector_2.
  Vector_2 base;

public:

  typedef FT Cartesian_coordinate_type;
  typedef const RT& Homogeneous_coordinate_type;
  typedef typename Vector_2::Cartesian_const_iterator Cartesian_const_iterator;
  typedef R_                                    R;

    PointH2() {}

    PointH2(const Origin &)
      : base(NULL_VECTOR) {}

    template < typename Tx, typename Ty >
    PointH2(const Tx & x, const Ty & y,
            typename boost::enable_if< boost::mpl::and_<boost::is_convertible<Tx, RT>,
                                                        boost::is_convertible<Ty, RT> > >::type* = 0)
      : base(x, y) {}

    PointH2(const FT& x, const FT& y)
      : base(x, y) {}

    PointH2(const RT& hx, const RT& hy, const RT& hw)
      : base(hx, hy, hw) {}

    bool    operator==( const PointH2<R>& p) const;
    bool    operator!=( const PointH2<R>& p) const;

    const RT & hx() const { return base.hx(); }
    const RT & hy() const { return base.hy(); }
    const RT & hw() const { return base.hw(); }

    FT      x()  const { return FT(hx()) / FT(hw()); }
    FT      y()  const { return FT(hy()) / FT(hw()); }

    FT      cartesian(int i)   const;
    FT      operator[](int i)  const;
    const RT & homogeneous(int i) const;

    Cartesian_const_iterator cartesian_begin() const
    {
      return base.cartesian_begin();
    }

    Cartesian_const_iterator cartesian_end() const
    {
      return base.cartesian_end();
    }

    int     dimension() const;

    Direction_2 direction() const;
};

template < class R >
inline
bool
PointH2<R>::operator==( const PointH2<R>& p) const
{
  return base == p.base;
}

template < class R >
inline
bool
PointH2<R>::operator!=( const PointH2<R>& p) const
{ return !(*this == p); }

template < class R >
inline
typename PointH2<R>::FT
PointH2<R>::cartesian(int i) const
{
  return base.cartesian(i);
}

template < class R >
inline
const typename PointH2<R>::RT &
PointH2<R>::homogeneous(int i) const
{
  return base.homogeneous(i);
}

template < class R >
inline
typename PointH2<R>::FT
PointH2<R>::operator[](int i) const
{ return base[i]; }


template < class R >
inline
int
PointH2<R>::dimension() const
{ return base.dimension(); }

template < class R >
inline
typename PointH2<R>::Direction_2
PointH2<R>::direction() const
{ return typename PointH2<R>::Direction_2(*this); }

} //namespace CGAL

#endif // CGAL_HOMOGENEOUS_POINT_2_H
