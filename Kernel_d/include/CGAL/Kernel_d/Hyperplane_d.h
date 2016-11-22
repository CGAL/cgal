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

#ifndef CGAL_HYPERPLANE_D_H
#define CGAL_HYPERPLANE_D_H

#include <CGAL/Dimension.h>
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Kernel_d/Direction_d.h>
#include <CGAL/enum.h>

namespace CGAL {

template <class pR>
class Hyperplane_d : public pR::Hyperplane_d_base
{
public:

  typedef CGAL::Dynamic_dimension_tag Ambient_dimension;
  typedef CGAL::Dynamic_dimension_tag Feature_dimension;

  typedef typename pR::Hyperplane_d_base Base;
  typedef Hyperplane_d<pR>               Self;
  typedef pR R;
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  typedef typename R::LA LA;

  Hyperplane_d(int d=0) : Base(d) {}
  Hyperplane_d(int a, int b, int c) : Base(a,b,c) {}
  Hyperplane_d(const RT& a, const RT& b, const RT& c) : 
    Base(a,b,c) {}
  Hyperplane_d(int a, int b, int c, int d) : Base(a,b,c,d) {}
  Hyperplane_d(const RT& a, const RT& b, const RT& c, const RT& d) : 
    Base(a,b,c,d) {}

  Hyperplane_d(const Point_d<R>& p, const Direction_d<R>& dir) :
    Base(p,dir) {}

  Hyperplane_d(const Hyperplane_d<R> &h) : Base(h) {}
  Hyperplane_d(const Base& p) : Base(p) {}

  template <class InputIterator>
  Hyperplane_d(int d, InputIterator first, InputIterator last)
    : Base (d, first, last) {}

  template <class InputIterator>
  Hyperplane_d(int d, InputIterator first, InputIterator last,
               const RT& D)
    : Base (d, first, last, D) {}

  template <class ForwardIterator>
  Hyperplane_d(ForwardIterator first, ForwardIterator last, 
               const Point_d<R>& o, Oriented_side side = ON_ORIENTED_BOUNDARY)
    : Base(first,last,o,side) {}

  Vector_d<R> orthogonal_vector() const 
  { return Base::orthogonal_vector(); }
  Direction_d<R> orthogonal_direction() const 
  { return Base::orthogonal_direction(); }

  bool operator==(const Self& w) const
  { return Base::operator==(w); }
  bool operator!=(const Self& w) const
  { return Base::operator!=(w); }
};

} //namespace CGAL

#endif //CGAL_HYPERPLANE_D_H
