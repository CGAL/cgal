// Copyright (c) 2001-2004  ENS of Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_CIRCLE_BY_RADIUS_2_H
#define CGAL_VISIBILITY_COMPLEX_2_CIRCLE_BY_RADIUS_2_H

#include <CGAL/Circle_2.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Visibility_complex_2/Rounded_sqrt.h>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {
// -----------------------------------------------------------------------------


template <class R_>
class Circle_by_radius_2 : public R_::Circle_2
{
public:
    // -------------------------------------------------------------------------
    typedef  R_                   R;
    typedef typename R::RT        RT;
    typedef typename R::FT        FT;
    typedef typename R::Circle_2  Circle_2;
    typedef typename R::Point_2   Point_2;
    // -------------------------------------------------------------------------
private:
    typename Rounded_sqrt<FT>::Sqrt sqrt;
public:
    Circle_by_radius_2()
	: Circle_2() , radius_(0) { }
    Circle_by_radius_2(const Circle_2 &t)
	: Circle_2(t) , radius_(sqrt(t.squared_radius())) { }
    Circle_by_radius_2(const Point_2&center,
		       const FT &radius,
		       const Orientation &orientation)
      : Circle_2(center, radius*radius, orientation) , radius_(radius) { }
    Circle_by_radius_2(const Point_2&center,
		       const FT &radius)
      : Circle_2(center, radius*radius, COUNTERCLOCKWISE) , radius_(radius) { }
    Circle_by_radius_2(const Point_2&p, const Point_2&q, const Point_2&r)
      : Circle_2(p,q,r) , 
	radius_(sqrt(Circle_2::squared_radius())) { }
    Circle_by_radius_2(const Point_2& p, const Point_2& q,
		       const Orientation &orientation)
      : Circle_2(p,q,orientation) , 
	radius_(sqrt(Circle_2::squared_radius())) { }
    Circle_by_radius_2(const Point_2& p,
		       const Point_2& q)
      : Circle_2(p,q,COUNTERCLOCKWISE) ,
	radius_(sqrt(Circle_2::squared_radius())) { }
    Circle_by_radius_2(const Point_2& center,
		       const Orientation& orientation)
      : Circle_2(center,FT(0),orientation) ,
	radius_(sqrt(Circle_2::squared_radius())) { }
    Circle_by_radius_2(const Point_2& center)
      : Circle_2(center,FT(0),COUNTERCLOCKWISE) ,
	radius_(sqrt(Circle_2::squared_radius())) { }
    // -------------------------------------------------------------------------
    FT radius() const { return radius_; }
private:
    FT radius_;
};

template < class R >
std::istream &
operator>>(std::istream &is, Circle_by_radius_2<R> &c)
{
  Circle_2<R> cc;
  is>>cc;
  if (is) c=Circle_by_radius_2<R>(cc);
  return is;
}


// -----------------------------------------------------------------------------
}
CGAL_END_NAMESPACE

#endif  // CGAL_VISIBILITY_COMPLEX_2_CIRCLE_BY_RADIUS_2_H
