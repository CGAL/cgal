// Copyright (c) 1998,2003  
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
// Author(s)     : Geert-Jan Giezeman, Sylvain Pion

#ifndef CGAL_KERNEL_WUTILS_H
#define CGAL_KERNEL_WUTILS_H

#include <CGAL/representation_tags.h>
#include <CGAL/determinant.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>

// Functions wmult() and wcross().

namespace CGAL {

namespace internal {

template < typename Rep_Tag > struct wmult_tag;

template <>
struct wmult_tag<Cartesian_tag>
{
  template < typename RT >
  const RT & operator()(const RT &a, const RT &) const
  { return a; }

  template < typename RT >
  const RT & operator()(const RT &a, const RT &, const RT &) const
  { return a; }

  template < typename RT >
  const RT & operator()(const RT &a, const RT &, const RT &, const RT &) const
  { return a; }

  template < typename RT >
  const RT & operator()(const RT &a, const RT &, const RT &, const RT &,
                        const RT &) const
  { return a; }
};

template <>
struct wmult_tag<Homogeneous_tag>
{
  template < typename RT >
  RT operator()(const RT &a, const RT &w) const
  { return a*w; }

  template < typename RT >
  RT operator()(const RT &a, const RT &w1, const RT &w2) const
  { return a*w1*w2; }

  template < typename RT >
  RT operator()(const RT &a, const RT &w1, const RT &w2, const RT &w3) const
  { return a*w1*w2*w3; }

  template < typename RT >
  RT operator()(const RT &a, const RT &w1, const RT &w2, const RT &w3,
                const RT &w4) const
  { return a*w1*w2*w3*w4; }
};

template < typename K >
struct wmult_functor
  : public wmult_tag<typename K::Rep_tag> {};


template < typename Rep_Tag > struct wmult_hw_tag;

template <>
struct wmult_hw_tag<Cartesian_tag>
{
  template < typename RT, typename T >
  const RT & operator()(const RT &a, const T &) const
  { return a; }

  template < typename RT, typename T >
  const RT & operator()(const RT &a, const RT &, const T &) const
  { return a; }

  template < typename RT, typename T >
  const RT & operator()(const RT &a, const RT &, const RT &, const T &) const
  { return a; }

  template < typename RT, typename T >
  const RT & operator()(const RT &a, const RT &, const RT &, const RT &,
                        const T &) const
  { return a; }
};

template <>
struct wmult_hw_tag<Homogeneous_tag>
{
  template < typename RT, typename T >
  RT operator()(const RT &a, const T &t) const
  { return a*t.hw(); }

  template < typename RT, typename T >
  RT operator()(const RT &a, const RT &w1, const T &t) const
  { return a*w1*t.hw(); }

  template < typename RT, typename T >
  RT operator()(const RT &a, const RT &w1, const RT &w2, const T &t) const
  { return a*w1*w2*t.hw(); }

  template < typename RT, typename T >
  RT operator()(const RT &a, const RT &w1, const RT &w2, const RT &w3,
                const T &t) const
  { return a*w1*w2*w3*t.hw(); }
};


template < typename K >
struct wmult_hw_functor
  : public wmult_hw_tag<typename K::Rep_tag> {};


template < typename Rep_Tag > struct wcross_tag_2;

template <>
struct wcross_tag_2<Cartesian_tag>
{
  template < typename Point_2 >
  typename Point_2::R::RT operator()(const Point_2 &p,
                                     const Point_2 &q,
                                     const Point_2 &r) const
  {
    return (q.x()-p.x())*(r.y()-q.y()) - (q.y()-p.y())*(r.x()-q.x());
  }
};

template <>
struct wcross_tag_2<Homogeneous_tag>
{
  template < typename Point_2 >
  typename Point_2::R::RT operator()(const Point_2 &p,
                                     const Point_2 &q,
                                     const Point_2 &r) const
  {
    return determinant(p.hx(), q.hx(), r.hx(),
                             p.hy(), q.hy(), r.hy(),
                             p.hw(), q.hw(), r.hw());
  }
};

template < typename K >
struct wcross_functor_2
  : public wcross_tag_2<typename K::Rep_tag> {};


template < typename Rep_Tag > struct wcross_tag_3;

template <>
struct wcross_tag_3<Cartesian_tag>
{
  template < typename Point_3 >
  typename Point_3::R::Vector_3 operator()(const Point_3 &p,
                                           const Point_3 &q,
                                           const Point_3 &r) const
  {
    typedef typename Point_3::R::FT         FT;
    typedef typename Point_3::R::Vector_3   Vector_3;
    FT x = (q.y()-p.y())*(r.z()-q.z()) - (q.z()-p.z())*(r.y()-q.y());
    FT y = (q.z()-p.z())*(r.x()-q.x()) - (q.x()-p.x())*(r.z()-q.z());
    FT z = (q.x()-p.x())*(r.y()-q.y()) - (q.y()-p.y())*(r.x()-q.x());
    return Vector_3(x, y, z);
  }
};

template <>
struct wcross_tag_3<Homogeneous_tag>
{
  template < typename Point_3 >
  typename Point_3::R::Vector_3 operator()(const Point_3 &p,
                                           const Point_3 &q,
                                           const Point_3 &r) const
  {
    typedef typename Point_3::R::RT         RT;
    typedef typename Point_3::R::Vector_3   Vector_3;
    RT x =  p.hy() * (q.hz()*r.hw() - q.hw()*r.hz() )
          + p.hz() * (q.hw()*r.hy() - q.hy()*r.hw() )
          + p.hw() * (q.hy()*r.hz() - q.hz()*r.hy() );
    RT y =  p.hz() * (q.hx()*r.hw() - q.hw()*r.hx() )
          + p.hx() * (q.hw()*r.hz() - q.hz()*r.hw() )
          + p.hw() * (q.hz()*r.hx() - q.hx()*r.hz() );
    RT z =  p.hx() * (q.hy()*r.hw() - q.hw()*r.hy() )
          + p.hy() * (q.hw()*r.hx() - q.hx()*r.hw() )
          + p.hw() * (q.hx()*r.hy() - q.hy()*r.hx() );
    return Vector_3(x, y, z);
  }
};

template < typename K >
struct wcross_functor_3
  : public wcross_tag_3<typename K::Rep_tag> {};

} // end namespace internal


// wmult_hw() is like wmult(), except it calls .hw() on its last argument.
// This way, we can completely avoid creating FT(1) for Cartesian.

template < typename K, typename T >
inline
typename K::RT
wmult_hw(K*, const typename K::RT &a,
             const T &t)
{
    return internal::wmult_hw_functor<K>()(a, t);
}

template < typename K, typename T >
inline
typename K::RT
wmult_hw(K*, const typename K::RT &a,
             const typename K::RT &w1,
             const T &t)
{
    return internal::wmult_hw_functor<K>()(a, w1, t);
}

template < typename K, typename T >
inline
typename K::RT
wmult_hw(K*, const typename K::RT &a,
             const typename K::RT &w1,
             const typename K::RT &w2,
             const T &t)
{
    return internal::wmult_hw_functor<K>()(a, w1, w2, t);
}

template < typename K, typename T >
inline
typename K::RT
wmult_hw(K*, const typename K::RT &a,
             const typename K::RT &w1,
             const typename K::RT &w2,
             const typename K::RT &w3,
             const T &t)
{
    return internal::wmult_hw_functor<K>()(a, w1, w2, w3, t);
}


template < typename K >
inline
typename K::RT
wmult(K*, const typename K::RT &a,
          const typename K::RT &w)
{
    return internal::wmult_functor<K>()(a, w);
}

template < typename K >
inline
typename K::RT
wmult(K*, const typename K::RT &a,
          const typename K::RT &w1,
          const typename K::RT &w2)
{
    return internal::wmult_functor<K>()(a, w1, w2);
}

template < typename K >
inline
typename K::RT
wmult(K*, const typename K::RT &a,
          const typename K::RT &w1,
          const typename K::RT &w2,
          const typename K::RT &w3)
{
    return internal::wmult_functor<K>()(a, w1, w2, w3);
}

template < typename K >
inline
typename K::RT
wmult(K*, const typename K::RT &a,
          const typename K::RT &w1,
          const typename K::RT &w2,
          const typename K::RT &w3,
          const typename K::RT &w4)
{
    return internal::wmult_functor<K>()(a, w1, w2, w3, w4);
}

template < typename K >
inline
typename K::RT
wcross(K*, const Point_2<K> &p,
           const Point_2<K> &q,
           const Point_2<K> &r)
{
    return internal::wcross_functor_2<K>()(p, q, r);
}

template < typename K >
inline
typename K::Vector_3
wcross(const Point_3<K> &p,
       const Point_3<K> &q,
       const Point_3<K> &r)
{
    return internal::wcross_functor_3<K>()(p, q, r);
}

} //namespace CGAL

#endif // CGAL_KERNEL_WUTILS_H
