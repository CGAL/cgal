// Copyright (c) 2000  
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
// Author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_ISO_CUBOID_3_H
#define CGAL_CARTESIAN_ISO_CUBOID_3_H

#include <CGAL/array.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Cartesian/predicates_on_points_3.h>

namespace CGAL {

template < class R_ >
class Iso_cuboidC3
{
  typedef typename R_::FT                   FT;
  typedef typename R_::Iso_cuboid_3         Iso_cuboid_3;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;
  typedef typename R_::Construct_point_3    Construct_point_3;

  typedef cpp11::array<Point_3, 2>          Rep;
  typedef typename R_::template Handle<Rep>::type  Base;

  Base base;

public:
  typedef R_                                R;

  Iso_cuboidC3() {}

  Iso_cuboidC3(const Point_3 &p, const Point_3 &q, int)
    : base(CGAL::make_array(p, q))
  {
    // I have to remove the assertions, because of Cartesian_converter.
    // CGAL_kernel_assertion(p.x()<=q.x());
    // CGAL_kernel_assertion(p.y()<=q.y());
    // CGAL_kernel_assertion(p.z()<=q.z());
  }

  Iso_cuboidC3(const Point_3 &p, const Point_3 &q)
  {
    Construct_point_3 construct_point_3;
    FT minx, maxx, miny, maxy, minz, maxz;
    if (p.x() < q.x()) { minx = p.x(); maxx = q.x(); }
    else               { minx = q.x(); maxx = p.x(); }
    if (p.y() < q.y()) { miny = p.y(); maxy = q.y(); }
    else               { miny = q.y(); maxy = p.y(); }
    if (p.z() < q.z()) { minz = p.z(); maxz = q.z(); }
    else               { minz = q.z(); maxz = p.z(); }
    base = Rep(CGAL::make_array(construct_point_3(minx, miny, minz),
	                         construct_point_3(maxx, maxy, maxz)));
  }

  Iso_cuboidC3(const Point_3 &left,   const Point_3 &right,
               const Point_3 &bottom, const Point_3 &top,
               const Point_3 &far_,   const Point_3 &close)
    : base(CGAL::make_array(Construct_point_3()(left.x(),  bottom.y(), far_.z()),
                             Construct_point_3()(right.x(), top.y(),    close.z())))
  {
    CGAL_kernel_precondition(!less_x(right, left));
    CGAL_kernel_precondition(!less_y(top, bottom));
    CGAL_kernel_precondition(!less_z(close, far_));
  }

  Iso_cuboidC3(const FT& min_x, const FT& min_y, const FT& min_z,
               const FT& max_x, const FT& max_y, const FT& max_z)
    : base(CGAL::make_array(Construct_point_3()(min_x, min_y, min_z),
	                     Construct_point_3()(max_x, max_y, max_z)))
  {
    CGAL_kernel_precondition(min_x <= max_x);
    CGAL_kernel_precondition(min_y <= max_y);
    CGAL_kernel_precondition(min_z <= max_z);
  }

  Iso_cuboidC3(const FT& min_hx, const FT& min_hy, const FT& min_hz,
               const FT& max_hx, const FT& max_hy, const FT& max_hz, 
               const FT& hw)
  {
    if (hw == FT(1))
       base = Rep(CGAL::make_array(Construct_point_3()(min_hx, min_hy, min_hz),
		                    Construct_point_3()(max_hx, max_hy, max_hz)));
    else
       base = Rep(CGAL::make_array(Construct_point_3()(min_hx/hw, min_hy/hw, min_hz/hw),
                                    Construct_point_3()(max_hx/hw, max_hy/hw, max_hz/hw)));
  }

  typename R::Boolean   operator==(const Iso_cuboidC3& s) const;
  typename R::Boolean   operator!=(const Iso_cuboidC3& s) const;

  const Point_3 & min BOOST_PREVENT_MACRO_SUBSTITUTION () const
  {
      return get_pointee_or_identity(base)[0];
  }
  const Point_3 & max BOOST_PREVENT_MACRO_SUBSTITUTION () const
  {
      return get_pointee_or_identity(base)[1];
  }
  Point_3 vertex(int i) const;
  Point_3 operator[](int i) const;

  Iso_cuboid_3 transform(const Aff_transformation_3 &t) const
  {
    return Iso_cuboidC3(t.transform((this->min)()), t.transform((this->max)()));
  }

  Bounded_side bounded_side(const Point_3& p) const;
  typename R::Boolean           has_on(const Point_3& p) const;
  typename R::Boolean           has_on_boundary(const Point_3& p) const;
  typename R::Boolean           has_on_bounded_side(const Point_3& p) const;
  typename R::Boolean           has_on_unbounded_side(const Point_3& p) const;
  typename R::Boolean           is_degenerate() const;
  const FT &   xmin() const;
  const FT &   ymin() const;
  const FT &   zmin() const;
  const FT &   xmax() const;
  const FT &   ymax() const;
  const FT &   zmax() const;
  const FT &   min_coord(int i) const;
  const FT &   max_coord(int i) const;

  FT           volume() const;
};

template < class R >
CGAL_KERNEL_INLINE
typename R::Boolean
Iso_cuboidC3<R>::operator==(const Iso_cuboidC3<R>& r) const
{ // FIXME : predicate
  if (CGAL::identical(base, r.base))
      return true;
  return (this->min)() == (r.min)() && (this->max)() == (r.max)();
}

template < class R >
inline
typename R::Boolean
Iso_cuboidC3<R>::operator!=(const Iso_cuboidC3<R>& r) const
{
  return !(*this == r);
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::xmin() const
{
  return (this->min)().x();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::ymin() const
{
  return (this->min)().y();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::zmin() const
{
  return (this->min)().z();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::xmax() const
{
  return (this->max)().x();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::ymax() const
{
  return (this->max)().y();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::zmax() const
{
  return (this->max)().z();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::min_coord(int i) const
{
  CGAL_kernel_precondition( i == 0 || i == 1 || i == 2 );
  if (i == 0)
     return xmin();
  else if (i == 1)
     return ymin();
  else 
     return zmin();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::max_coord(int i) const
{
  CGAL_kernel_precondition( i == 0 || i == 1 || i == 2 );
  if (i == 0)
     return xmax();
  else if (i == 1)
     return ymax();
  else 
     return zmax();
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Iso_cuboidC3<R>::Point_3
Iso_cuboidC3<R>::vertex(int i) const
{
  Construct_point_3 construct_point_3;
  switch (i%8)
  {
    case 0: return (this->min)();
    case 1: return construct_point_3((this->max)().hx(), (this->min)().hy(), (this->min)().hz());
    case 2: return construct_point_3((this->max)().hx(), (this->max)().hy(), (this->min)().hz());
    case 3: return construct_point_3((this->min)().hx(), (this->max)().hy(), (this->min)().hz());
    case 4: return construct_point_3((this->min)().hx(), (this->max)().hy(), (this->max)().hz());
    case 5: return construct_point_3((this->min)().hx(), (this->min)().hy(), (this->max)().hz());
    case 6: return construct_point_3((this->max)().hx(), (this->min)().hy(), (this->max)().hz());
    default: // case 7:
        return (this->max)();
  }
}

template < class R >
inline
typename Iso_cuboidC3<R>::Point_3
Iso_cuboidC3<R>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
inline
typename Iso_cuboidC3<R>::FT
Iso_cuboidC3<R>::volume() const
{
  return (xmax()-xmin()) * (ymax()-ymin()) * (zmax()-zmin());
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
Iso_cuboidC3<R>::
bounded_side(const typename Iso_cuboidC3<R>::Point_3& p) const
{
  if (strict_dominance(p, (this->min)()) && strict_dominance((this->max)(), p) )
    return ON_BOUNDED_SIDE;
  if (dominance(p, (this->min)()) && dominance((this->max)(), p))
    return ON_BOUNDARY;
  return ON_UNBOUNDED_SIDE;
}

template < class R >
inline
typename R::Boolean
Iso_cuboidC3<R>::
has_on_boundary(const typename Iso_cuboidC3<R>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
typename R::Boolean
Iso_cuboidC3<R>::
has_on(const typename Iso_cuboidC3<R>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
typename R::Boolean
Iso_cuboidC3<R>::
has_on_bounded_side(const typename Iso_cuboidC3<R>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
typename R::Boolean
Iso_cuboidC3<R>::
has_on_unbounded_side(const typename Iso_cuboidC3<R>::Point_3& p)
    const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
typename R::Boolean
Iso_cuboidC3<R>::is_degenerate() const
{ // FIXME : predicate
  return (this->min)().hx() == (this->max)().hx()
      || (this->min)().hy() == (this->max)().hy()
      || (this->min)().hz() == (this->max)().hz();
}

} //namespace CGAL

#endif // CGAL_CARTESIAN_ISO_CUBOID_3_H
