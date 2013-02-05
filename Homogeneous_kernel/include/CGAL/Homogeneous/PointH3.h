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
//
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_HOMOGENEOUS_POINT_3_H
#define CGAL_HOMOGENEOUS_POINT_3_H

#include <CGAL/Origin.h>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/logical.hpp>

namespace CGAL {

template < class R_ >
class PointH3
{
   typedef typename R_::RT                   RT;
   typedef typename R_::FT                   FT;
   typedef typename R_::Vector_3             Vector_3;
   typedef typename R_::Point_3              Point_3;
   typedef typename R_::Direction_3          Direction_3;
   typedef typename R_::Aff_transformation_3 Aff_transformation_3;

   typedef Rational_traits<FT>  Rat_traits;

   // Reference-counting is handled in Vector_3.
   Vector_3 base;

public:

   typedef typename Vector_3::Cartesian_const_iterator Cartesian_const_iterator;
   typedef R_                 R;

  PointH3() {}

  PointH3(const Origin &)
    : base(NULL_VECTOR) {}

  template < typename Tx, typename Ty, typename Tz >
  PointH3(const Tx & x, const Ty & y, const Tz & z,
          typename boost::enable_if< boost::mpl::and_< boost::mpl::and_< boost::is_convertible<Tx, RT>,
                                                                         boost::is_convertible<Ty, RT> >,
                                                       boost::is_convertible<Tz, RT> > >::type* = 0)
    : base(x, y, z) {}

  PointH3(const FT& x, const FT& y, const FT& z)
    : base(x, y, z) {}

  PointH3(const RT& x, const RT& y, const RT& z, const RT& w)
    : base(x, y, z, w) {}

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
    return base.cartesian_begin();
  }

  Cartesian_const_iterator cartesian_end() const
  {
    return base.cartesian_end();
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
{ return base.hx(); }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hy() const
{ return base.hy(); }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hz() const
{ return base.hz(); }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hw() const
{ return base.hw(); }

template < class R >
inline
typename PointH3<R>::FT
PointH3<R>::x()  const
{ return base.x(); }

template < class R >
inline
typename PointH3<R>::FT
PointH3<R>::y()  const
{ return base.y(); }

template < class R >
inline
typename PointH3<R>::FT
PointH3<R>::z()  const
{ return base.z(); }

template < class R >
inline
int
PointH3<R>::dimension() const
{ return base.dimension(); }

template < class R >
inline
typename PointH3<R>::FT
PointH3<R>::cartesian(int i) const
{
  return base.cartesian(i);
}

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::homogeneous(int i) const
{
  return base.homogeneous(i);
}

template < class R >
inline
typename PointH3<R>::FT
PointH3<R>::operator[](int i) const
{
  return base[i];
}

template < class R >
inline
typename PointH3<R>::Direction_3
PointH3<R>::direction() const
{ return Direction_3(*this); }

template < class R >
inline
bool
PointH3<R>::operator==( const PointH3<R> & p) const
{
  return base == p.base;
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

} //namespace CGAL

#endif // CGAL_HOMOGENEOUS_POINT_3_H
