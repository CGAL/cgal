// Copyright (c) 2000,2001
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
// Author(s)     : Michael Seel

#ifndef CGAL_POINT_D_H
#define CGAL_POINT_D_H

#include <CGAL/Dimension.h>
#include <CGAL/Origin.h>

namespace CGAL {


template <class pR>
  class Direction_d;

template <class pR>
  class Vector_d;


template <class pR>
class Point_d : public pR::Point_d_base
{ public:
  typedef typename pR::Point_d_base Base;
  typedef Point_d<pR>               Self;
  typedef pR R;
private:
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  typedef typename R::LA LA;
public:

  typedef CGAL::Dynamic_dimension_tag Ambient_dimension;
  typedef CGAL::Dimension_tag<0>      Feature_dimension;
    template < typename Kernel2 >
        struct WithAnotherKernel
        {
            typedef Point_d<Kernel2>  Type;
        };

  Point_d(int d=0) : Base(d) {}
  Point_d(int d, const Origin &o) : Base(d,o) {}

  Point_d(int a, int b, int c = 1) :
    Base(RT(a),RT(b),RT(c)) {}
  Point_d(const RT& a, const RT& b, const RT& c = 1) :
    Base(a,b,c) {}
  Point_d(int a, int b, int c, int d) :
    Base(RT(a),RT(b),RT(c),RT(d)) {}
  Point_d(const RT& a, const RT& b, const RT& c, const RT& d) :
    Base(a,b,c,d) {}

  template <class InputIterator>
  Point_d (int d, InputIterator first, InputIterator last)
    : Base (d, first, last) {}
  template <class InputIterator>
  Point_d(int d, InputIterator first, InputIterator last, const RT& D)
    : Base (d, first, last, D) {}

  Point_d(const Base& p) : Base(p) {}

  Vector_d<R> operator-(const Origin& o) const
  { return Base::operator-(o); }
  Vector_d<R> operator-(const Self& q) const
  { return Base::operator-(q); }
  Self operator+(const Vector_d<R>& v) const
  { return Base::operator+(v); }
  Self operator-(const Vector_d<R>& v) const
  { return Base::operator-(v); }
  Self& operator+=(const Vector_d<R>& v)
  { return static_cast<Self&>(Base::operator+=(v)); }
  Self& operator-=(const Vector_d<R>& v)
  { return static_cast<Self&>(Base::operator-=(v)); }

  inline bool operator<(const Self& q) const
  { return R().less_lexicographically_d_object()(*this, q); }
  inline bool operator>(const Self& q) const
  { return R().less_lexicographically_d_object()(q, *this); }
  inline bool operator<=(const Self& q) const
  { return ! R().less_lexicographically_d_object()(q, *this); }
  inline bool operator>=(const Self& q) const
  { return ! R().less_lexicographically_d_object()(*this, q); }
};

} //namespace CGAL
#endif //CGAL_POINT_D_H
