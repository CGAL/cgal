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

#ifndef CGAL_ISO_CUBOIDH3_H
#define CGAL_ISO_CUBOIDH3_H

#include <CGAL/array.h>
#include <CGAL/enum.h>
#include <CGAL/kernel_config.h>
#include <CGAL/kernel_assertions.h>

namespace CGAL {

template <class R_>
class Iso_cuboidH3
{
  typedef typename R_::RT                   RT;
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef cpp11::array<Point_3, 2>          Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                 R;

  Iso_cuboidH3() {}

  Iso_cuboidH3(const Point_3& p, const Point_3& q, int)
    : base(CGAL::make_array(p, q))
  {
    CGAL_kernel_assertion(p.x()<=q.x());
    CGAL_kernel_assertion(p.y()<=q.y());
    CGAL_kernel_assertion(p.z()<=q.z());
  }

  Iso_cuboidH3(const Point_3& p, const Point_3& q);

  Iso_cuboidH3(const Point_3& left,   const Point_3& right,
               const Point_3& bottom, const Point_3& top,
               const Point_3& far_,   const Point_3& close);

  Iso_cuboidH3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz, 
               const RT& hw);

  Iso_cuboidH3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz);

  bool      operator==(const Iso_cuboidH3<R>& s) const;
  bool      operator!=(const Iso_cuboidH3<R>& s) const;

  const Point_3 & min BOOST_PREVENT_MACRO_SUBSTITUTION () const;
  const Point_3 & max BOOST_PREVENT_MACRO_SUBSTITUTION () const;
  Point_3  vertex(int i) const;
  Point_3  operator[](int i) const;

  Iso_cuboidH3<R>
            transform(const Aff_transformation_3& t) const;
  Bounded_side
            bounded_side(const Point_3& p) const;
  bool      has_on(const Point_3& p) const;
  bool      has_on_boundary(const Point_3& p) const;
  bool      has_on_bounded_side(const Point_3& p) const;
  bool      has_on_unbounded_side(const Point_3& p) const;
  bool      is_degenerate() const;
  FT        xmin() const;
  FT        ymin() const;
  FT        zmin() const;
  FT        xmax() const;
  FT        ymax() const;
  FT        zmax() const;
  FT        min_coord(int i) const;
  FT        max_coord(int i) const;

  FT        volume() const;

};

template < class R >
CGAL_KERNEL_LARGE_INLINE
Iso_cuboidH3<R>::
Iso_cuboidH3(const typename Iso_cuboidH3<R>::Point_3& p,
	     const typename Iso_cuboidH3<R>::Point_3& q)
{
  bool px_k_qx = ( p.hx()*q.hw() < q.hx()*p.hw() );
  bool py_k_qy = ( p.hy()*q.hw() < q.hy()*p.hw() );
  bool pz_k_qz = ( p.hz()*q.hw() < q.hz()*p.hw() );

  RT minx;
  RT miny;
  RT minz;
  RT maxx;
  RT maxy;
  RT maxz;
  RT minw = p.hw()*q.hw();
  RT maxw = p.hw()*q.hw();
  if ( px_k_qx )
  {
      minx = p.hx()*q.hw();
      maxx = q.hx()*p.hw();
  }
  else
  {
      minx = q.hx()*p.hw();
      maxx = p.hx()*q.hw();
  }
  if ( py_k_qy )
  {
      miny = p.hy()*q.hw();
      maxy = q.hy()*p.hw();
  }
  else
  {
      miny = q.hy()*p.hw();
      maxy = p.hy()*q.hw();
  }
  if ( pz_k_qz )
  {
      minz = p.hz()*q.hw();
      maxz = q.hz()*p.hw();
  }
  else
  {
      minz = q.hz()*p.hw();
      maxz = p.hz()*q.hw();
  }
  base = Rep(CGAL::make_array(Point_3(minx, miny, minz, minw),
                               Point_3(maxx, maxy, maxz, maxw)));
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Iso_cuboidH3<R>::
Iso_cuboidH3(const typename Iso_cuboidH3<R>::Point_3& left,
             const typename Iso_cuboidH3<R>::Point_3& right,
             const typename Iso_cuboidH3<R>::Point_3& bottom,
             const typename Iso_cuboidH3<R>::Point_3& top,
             const typename Iso_cuboidH3<R>::Point_3& far_,
             const typename Iso_cuboidH3<R>::Point_3& close)
  : base(CGAL::make_array(Point_3(left.hx()   * bottom.hw() * far_.hw(),
                                   bottom.hy() * left.hw()   * far_.hw(),
                                   far_.hz()   * left.hw()   * bottom.hw(),
                                   left.hw()   * bottom.hw() * far_.hw()),
                           Point_3(right.hx()  * top.hw()    * close.hw(),
                                   top.hy()    * right.hw()  * close.hw(),
                                   close.hz()  * right.hw()  * top.hw(),
                                   right.hw()  * top.hw()    * close.hw())))
{
  CGAL_kernel_precondition(!less_x(right, left));
  CGAL_kernel_precondition(!less_y(top, bottom));
  CGAL_kernel_precondition(!less_z(close, far_));
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Iso_cuboidH3<R>::
Iso_cuboidH3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
             const RT& max_hx, const RT& max_hy, const RT& max_hz)
  : base(CGAL::make_array(Point_3(min_hx, min_hy, min_hz, RT(1)),
                           Point_3(max_hx, max_hy, max_hz, RT(1))))
{}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Iso_cuboidH3<R>::
Iso_cuboidH3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
             const RT& max_hx, const RT& max_hy, const RT& max_hz, 
             const RT& hw)
  : base(CGAL::make_array(Point_3(min_hx, min_hy, min_hz, hw),
                           Point_3(max_hx, max_hy, max_hz, hw)))
{}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidH3<R>::
operator==(const Iso_cuboidH3<R>& r) const
{ return  ((this->min)() == (r.min)()) && ((this->max)() == (r.max)()); }

template < class R >
inline
bool
Iso_cuboidH3<R>::
operator!=(const Iso_cuboidH3<R>& r) const
{ return !(*this == r); }

template < class R >
inline
const typename Iso_cuboidH3<R>::Point_3 &
Iso_cuboidH3<R>::min BOOST_PREVENT_MACRO_SUBSTITUTION () const
{ return get_pointee_or_identity(base)[0]; }

template < class R >
inline
const typename Iso_cuboidH3<R>::Point_3 &
Iso_cuboidH3<R>::max BOOST_PREVENT_MACRO_SUBSTITUTION () const
{ return get_pointee_or_identity(base)[1]; }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::xmin() const
{ return  FT( (this->min)().hx() ) / FT( (this->min)().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::ymin() const
{ return  FT( (this->min)().hy() ) / FT( (this->min)().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::zmin() const
{ return  FT( (this->min)().hz() ) / FT( (this->min)().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::xmax() const
{ return  FT( (this->max)().hx() ) / FT( (this->max)().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::ymax() const
{ return  FT( (this->max)().hy() ) / FT( (this->max)().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::zmax() const
{ return  FT( (this->max)().hz() ) / FT( (this->max)().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::min_coord(int i) const
{ 
   CGAL_kernel_precondition(i == 0 || i == 1 || i == 2);
   if ( i == 0 )
       return xmin();
   else if (i == 1)
       return ymin();
   return zmin();
}

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::max_coord(int i) const
{ 
   CGAL_kernel_precondition(i == 0 || i == 1 || i == 2);
   if ( i == 0 )
      return xmax();
   else if ( i == 1 )
      return ymax();
   return zmax();
}

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::volume() const
{ return  (xmax() - xmin()) * (ymax() - ymin()) * (zmax() - zmin()); }

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Iso_cuboidH3<R>::Point_3
Iso_cuboidH3<R>::vertex(int i) const
{
  switch (i%8)
  {
    case 0: return (this->min)();
    case 1: return Point_3( (this->max)().hx(), (this->min)().hy(),
		                   (this->min)().hz(), (this->min)().hw() );
    case 2: return Point_3( (this->max)().hx(), (this->max)().hy(),
		                   (this->min)().hz(), (this->min)().hw() );
    case 3: return Point_3( (this->min)().hx(), (this->max)().hy(),
		                   (this->min)().hz(), (this->min)().hw() );
    case 4: return Point_3( (this->min)().hx(), (this->max)().hy(),
		                   (this->max)().hz(), (this->min)().hw() );
    case 5: return Point_3( (this->min)().hx(), (this->min)().hy(),
		                   (this->max)().hz(), (this->min)().hw() );
    case 6: return Point_3( (this->max)().hx(), (this->min)().hy(),
		                   (this->max)().hz(), (this->min)().hw() );
    default: /*case 7:*/ return (this->max)();
  }
}

template < class R >
inline
typename Iso_cuboidH3<R>::Point_3
Iso_cuboidH3<R>::operator[](int i) const
{ return vertex(i); }

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
Iso_cuboidH3<R>::
bounded_side(const typename Iso_cuboidH3<R>::Point_3& p) const
{
  if (    (p.hx()*(this->min)().hw() < (this->min)().hx()*p.hw() )
        ||(p.hy()*(this->min)().hw() < (this->min)().hy()*p.hw() )
        ||(p.hz()*(this->min)().hw() < (this->min)().hz()*p.hw() )
        ||(p.hx()*(this->max)().hw() > (this->max)().hx()*p.hw() )
        ||(p.hy()*(this->max)().hw() > (this->max)().hy()*p.hw() )
        ||(p.hz()*(this->max)().hw() > (this->max)().hz()*p.hw() )  )
  { return ON_UNBOUNDED_SIDE; }
  if (    (p.hx()*(this->min)().hw() == (this->min)().hx()*p.hw() )
        ||(p.hy()*(this->min)().hw() == (this->min)().hy()*p.hw() )
        ||(p.hz()*(this->min)().hw() == (this->min)().hz()*p.hw() )
        ||(p.hx()*(this->max)().hw() == (this->max)().hx()*p.hw() )
        ||(p.hy()*(this->max)().hw() == (this->max)().hy()*p.hw() )
        ||(p.hz()*(this->max)().hw() == (this->max)().hz()*p.hw() )  )
  { return ON_BOUNDARY; }
  else
  { return ON_BOUNDED_SIDE; }
}

template < class R >
inline
bool
Iso_cuboidH3<R>::
has_on_boundary(const typename Iso_cuboidH3<R>::Point_3& p) const
{ return ( bounded_side(p) == ON_BOUNDARY ); }

template < class R >
inline
bool
Iso_cuboidH3<R>::has_on(const typename Iso_cuboidH3<R>::Point_3& p) const
{ return ( bounded_side(p) == ON_BOUNDARY ); }

template < class R >
inline
bool
Iso_cuboidH3<R>::
has_on_bounded_side(const typename Iso_cuboidH3<R>::Point_3& p) const
{ return ( bounded_side(p) == ON_BOUNDED_SIDE ); }

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidH3<R>::
has_on_unbounded_side(const typename Iso_cuboidH3<R>::Point_3& p) const
{
  return (   ( lexicographically_xyz_smaller(p,(this->min)() ))
           ||( lexicographically_xyz_smaller((this->max)(),p ))  );
}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidH3<R>::is_degenerate() const
{
  return (  ( (this->min)().hx() == (this->max)().hx() )
         || ( (this->min)().hy() == (this->max)().hy() )
         || ( (this->min)().hz() == (this->max)().hz() ) );
}

template < class R >
CGAL_KERNEL_INLINE
Iso_cuboidH3<R>
Iso_cuboidH3<R>::
transform(const typename Iso_cuboidH3<R>::Aff_transformation_3&t) const
{
  return Iso_cuboidH3<R>(t.transform((this->min)() ),
                             t.transform((this->max)() ) );
}

} //namespace CGAL

#endif // CGAL_ISO_CUBOIDH3_H
