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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_HOMOGENEOUS_POINT_3_H
#define CGAL_HOMOGENEOUS_POINT_3_H

#include <CGAL/Origin.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Fourtuple.h>
#include <CGAL/Kernel/Cartesian_coordinate_iterator_3.h>

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

   typedef Fourtuple<RT>                            Rep;
   typedef typename R_::template Handle<Rep>::type  Base;

   Base base;

public:
   typedef Cartesian_coordinate_iterator_3<R_> Cartesian_const_iterator;
   typedef R_                 R;

  PointH3() {}

  PointH3(const Origin &)
    : base (RT(0), RT(0), RT(0), RT(1)) { }

  PointH3(const RT& x, const RT& y, const RT& z)
    : base(x, y, z, RT(1)) {}

  PointH3(const RT& x, const RT& y, const RT& z, const RT& w)
  {
    if ( w < RT(0) )
      base = Rep(-x,-y,-z,-w);
    else
      base = Rep(x,y,z,w);
  }

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
  Bbox_3      bbox() const;

  bool  operator==( const PointH3<R>& p) const;
  bool  operator!=( const PointH3<R>& p) const;
};


template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hx() const
{ return get(base).e0 ; }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hy() const
{ return get(base).e1 ; }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hz() const
{ return get(base).e2 ; }

template < class R >
inline
const typename PointH3<R>::RT &
PointH3<R>::hw() const
{ return get(base).e3 ; }

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::FT
PointH3<R>::x()  const
{ return ( FT(hx()) / FT(hw())); }

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::FT
PointH3<R>::y()  const
{ return ( FT(hy()) / FT(hw())); }

template < class R >
CGAL_KERNEL_INLINE
typename PointH3<R>::FT
PointH3<R>::z()  const
{ return ( FT(hz()) / FT(hw())); }

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
  switch (i)
  {
      case 0:  return x();
      case 1:  return y();
  }
  return z();
}

template < class R >
CGAL_KERNEL_INLINE
const typename PointH3<R>::RT &
PointH3<R>::homogeneous(int i) const
{
  CGAL_kernel_precondition(i == 0 || i == 1 || i == 2 || i == 3);
  switch (i)
  {
     case 0:   return hx();
     case 1:   return hy();
     case 2:   return hz();
  }
  return hw();
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

#ifndef CGAL_NO_OSTREAM_INSERT_POINTH3
template < class R >
std::ostream &operator<<(std::ostream &os, const PointH3<R> &p)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << p.hx() << ' ' << p.hy() << ' ' << p.hz() << ' ' << p.hw();
    case IO::BINARY :
        write(os, p.hx());
        write(os, p.hy());
        write(os, p.hz());
        write(os, p.hw());
        return os;
    default:
        return os << "PointH3(" << p.hx() << ", "
                                << p.hy() << ", "
                                << p.hz() << ", "
                                << p.hw() << ')';
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_POINTH3

#ifndef CGAL_NO_ISTREAM_EXTRACT_POINTH3
template < class R >
std::istream &operator>>(std::istream &is, PointH3<R> &p)
{
  typename R::RT hx, hy, hz, hw;
  switch(is.iword(IO::mode)) {
  case IO::ASCII :
        is >> hx >> hy >> hz >> hw;
        break;
  case IO::BINARY :
        read(is, hx);
        read(is, hy);
        read(is, hz);
        read(is, hw);
        break;
  default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  p = PointH3<R>(hx, hy, hz, hw);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_POINTH3


template < class R >
inline
typename R::Point_3
PointH3<R>::transform(const typename PointH3<R>::Aff_transformation_3& t) const
{ return t.transform(static_cast<const Point_3&>(*this)); }

template < class R >
CGAL_KERNEL_LARGE_INLINE
Bbox_3
PointH3<R>::bbox() const
{
   Interval_nt<> ihx = CGAL_NTS to_interval(hx());
   Interval_nt<> ihy = CGAL_NTS to_interval(hy());
   Interval_nt<> ihz = CGAL_NTS to_interval(hz());
   Interval_nt<> ihw = CGAL_NTS to_interval(hw());

   Interval_nt<> ix = ihx/ihw;
   Interval_nt<> iy = ihy/ihw;
   Interval_nt<> iz = ihz/ihw;

   return Bbox_3(ix.inf(), iy.inf(), iz.inf(), ix.sup(), iy.sup(), iz.sup());
}

CGAL_END_NAMESPACE

#endif // CGAL_HOMOGENEOUS_POINT_3_H
