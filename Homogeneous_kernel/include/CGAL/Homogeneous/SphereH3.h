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

#ifndef CGAL_SPHEREH3_H
#define CGAL_SPHEREH3_H

#include <CGAL/Interval_nt.h>
#include <CGAL/Homogeneous/predicates_on_pointsH3.h>
#include <boost/tuple/tuple.hpp>
#include <CGAL/Kernel/global_functions_3.h>

namespace CGAL {

template <class R_>
class SphereH3
{
   typedef typename R_::RT                   RT;
   typedef typename R_::FT                   FT;
   typedef typename R_::Point_3              Point_3;

   typedef boost::tuple<Point_3, FT, Orientation>   Rep;
   typedef typename R_::template Handle<Rep>::type  Base;

   Base base;

public:
   typedef R_                R;

      SphereH3() {}

      SphereH3(const Point_3& p, const FT& sq_rad,
               const Orientation& o = COUNTERCLOCKWISE);

      SphereH3(const Point_3& p, const Point_3& q,
               const Point_3& r, const Point_3& u);

      SphereH3(const Point_3& p, const Point_3& q,
               const Point_3& r,
               const Orientation& o = COUNTERCLOCKWISE);

      SphereH3(const Point_3&  p, const Point_3&  q,
               const Orientation& o = COUNTERCLOCKWISE);

      SphereH3(const Point_3&  p,
               const Orientation& o = COUNTERCLOCKWISE);

      bool
      operator==(const SphereH3<R>&) const;

      bool
      operator!=(const SphereH3<R>& s) const
      { return !(*this == s); }

      const Point_3 & center() const;

      const FT & squared_radius() const;

      Orientation orientation() const;

      bool is_degenerate() const;

      SphereH3<R> opposite() const;

      Oriented_side oriented_side(const Point_3& p) const;

      bool
      has_on_boundary(const Point_3& p) const
      { return oriented_side(p)==ON_ORIENTED_BOUNDARY; }

      bool
      has_on_positive_side(const Point_3& p) const
      { return oriented_side(p)==ON_POSITIVE_SIDE; }

      bool
      has_on_negative_side(const Point_3& p) const
      { return oriented_side(p)==ON_NEGATIVE_SIDE; }

      Bounded_side
      bounded_side(const Point_3& p) const;

      bool
      has_on_bounded_side(const Point_3& p) const
      { return bounded_side(p)==ON_BOUNDED_SIDE; }

      bool
      has_on_unbounded_side(const Point_3& p) const
      { return bounded_side(p)==ON_UNBOUNDED_SIDE; }
};


template < class R >
CGAL_KERNEL_INLINE
SphereH3<R>::SphereH3(const typename SphereH3<R>::Point_3& center,
                      const FT& squared_radius,
                      const Orientation& o)
{
  CGAL_kernel_precondition( !( squared_radius < FT(0))
                          &&( o != COLLINEAR) );
  base = Rep(center, squared_radius, o);
}

template <class R>
CGAL_KERNEL_INLINE
SphereH3<R>::SphereH3(const typename SphereH3<R>::Point_3& center,
                      const Orientation& o)
{
  CGAL_kernel_precondition( ( o != COLLINEAR) );
  base = Rep(center, FT(0), o);
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
SphereH3<R>::SphereH3(const typename SphereH3<R>::Point_3& p,
                      const typename SphereH3<R>::Point_3& q,
                      const Orientation& o)
{
  CGAL_kernel_precondition( o != COLLINEAR);
  Point_3 center = midpoint(p,q);
  FT     squared_radius = squared_distance(p,center);
  base = Rep(center, squared_radius, o);
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
SphereH3<R>::SphereH3(const typename SphereH3<R>::Point_3& p,
                      const typename SphereH3<R>::Point_3& q,
                      const typename SphereH3<R>::Point_3& r,
                      const Orientation& o)
{
  CGAL_kernel_precondition( o != COLLINEAR);
  Point_3 center = circumcenter(p,q,r);
  FT     squared_radius = squared_distance(p,center);
  base = Rep(center, squared_radius, o);
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
SphereH3<R>::SphereH3(const typename SphereH3<R>::Point_3& p,
                      const typename SphereH3<R>::Point_3& q,
                      const typename SphereH3<R>::Point_3& r,
                      const typename SphereH3<R>::Point_3& s)
{
  Orientation o = CGAL::orientation(p,q,r,s);
  CGAL_kernel_precondition( o != COLLINEAR);
  Point_3 center = circumcenter(p,q,r,s);
  FT     squared_radius = squared_distance(p,center);
  base = Rep(center, squared_radius, o);
}

template <class R>
CGAL_KERNEL_INLINE
bool
SphereH3<R>::operator==(const SphereH3<R>& s) const
{
   return    ( orientation() == s.orientation())
          && ( center() == s.center())
          && ( squared_radius() == s.squared_radius());
}

template <class R>
inline
const typename SphereH3<R>::Point_3 &
SphereH3<R>::center() const
{ return get_pointee_or_identity(base).template get<0>(); }

template <class R>
inline
const typename SphereH3<R>::FT &
SphereH3<R>::squared_radius() const
{ return get_pointee_or_identity(base).template get<1>(); }

template <class R>
inline
Orientation
SphereH3<R>::orientation() const
{ return get_pointee_or_identity(base).template get<2>(); }

template <class R>
inline
bool
SphereH3<R>::is_degenerate() const
{ return squared_radius() <= FT(0) ; }

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
SphereH3<R>::oriented_side(const typename SphereH3<R>::Point_3& p) const
{ return Oriented_side(bounded_side(p) * orientation()); }

template <class R>
CGAL_KERNEL_INLINE
Bounded_side
SphereH3<R>::bounded_side(const typename SphereH3<R>::Point_3& p) const
{
  return Bounded_side(CGAL_NTS compare(squared_radius(),
                                       squared_distance(center(),p)));
}

template <class R>
inline
SphereH3<R>
SphereH3<R>::opposite() const
{
  return SphereH3<R>(center(), squared_radius(),
                         CGAL::opposite(orientation()) );
}

} //namespace CGAL

#endif // CGAL_SPHEREH3_H
