// Copyright (c) 2000  Utrecht University (The Netherlands),
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
// $URL$
// $Id$
// 
//
// Author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_POINT_2_H
#define CGAL_CARTESIAN_POINT_2_H

#include <CGAL/Origin.h>

namespace CGAL {

template < class R_ >
class PointC2
{
  typedef PointC2<R_>                       Self;
  typedef typename R_::FT                   FT;
  typedef typename R_::Vector_2             Vector_2;
  typedef typename R_::Point_2              Point_2;

  // We do not use reference counting here as it is done at the Vector_2 level.
  Vector_2 base;

public:

  typedef typename Vector_2::Cartesian_const_iterator Cartesian_const_iterator;
  
  typedef R_                                R;

  PointC2() {}

  PointC2(const Origin &)
    : base(NULL_VECTOR) {}

  PointC2(const FT &x, const FT &y)
    : base(x, y) {}

  PointC2(const FT &hx, const FT &hy, const FT &hw)
    : base(hx, hy, hw) {}

  const FT& x() const
  {
      return base.x();
  }
  
  const FT& y() const
  {
      return base.y();
  }

  const FT& hx() const
  {
      return base.hx();
  }
  const FT& hy() const
  {
      return base.hy();
  }
  const FT& hw() const
  {
      return base.hw();
  }

  Cartesian_const_iterator cartesian_begin() const 
  {
    return base.cartesian_begin(); 
  }

  Cartesian_const_iterator cartesian_end() const 
  {
    return base.cartesian_end(); 
  }

  typename R_::Boolean   operator==(const PointC2 &p) const
  {
      return base == p.base;
  }
  typename R_::Boolean   operator!=(const PointC2 &p) const
  {
      return !(*this == p);
  }

};

} //namespace CGAL

#endif // CGAL_CARTESIAN_POINT_2_H
