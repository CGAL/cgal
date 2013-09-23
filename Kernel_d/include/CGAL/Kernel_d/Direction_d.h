// Copyright (c) 2000,2001  
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
// Author(s)     : Michael Seel

#ifndef CGAL_DIRECTION_D_H
#define CGAL_DIRECTION_D_H

#include <CGAL/Dimension.h>

namespace CGAL {

template <class pR>
  class Vector_d;

template <class pR>
class Direction_d : public pR::Direction_d_base
{ public:

  typedef CGAL::Dynamic_dimension_tag Ambient_dimension;
  typedef CGAL::Dimension_tag<0>      Feature_dimension;

  typedef typename pR::Direction_d_base Base;
  typedef Direction_d<pR>               Self;
  typedef pR R;
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  typedef typename R::LA LA;
  typedef typename Base::Base_direction Base_direction;

  Direction_d(int d=0) : Base(d) {}
  Direction_d(int a, int b) : Base(a,b) {}
  Direction_d(const RT& a, const RT& b) : Base(a,b) {}
  Direction_d(int a, int b, int c) : Base(a,b,c) {}
  Direction_d(const RT& a, const RT& b, const RT& c) : Base(a,b,c) {}

  template <class InputIterator>
  Direction_d (int d, InputIterator first, InputIterator last)
    : Base(d, first, last) {}

  Direction_d(const Direction_d<R> &d) : Base(d) {}
  Direction_d(const Vector_d<R> &v) : Base(v) {}
  Direction_d(int d, Base_direction, int i) : 
    Base(d,Base_direction(),i) {}
  Direction_d(const Base& p) : Base(p) {}

  Self operator-() const { return Base::operator-(); }
  Vector_d<R> vector() const { return Base::vector(); }

  bool operator==(const Self& w) const
  { return Base::operator==(w); }
  bool operator!=(const Self& w) const
  { return Base::operator!=(w); }
};

} //namespace CGAL
#endif //CGAL_DIRECTION_D_H
