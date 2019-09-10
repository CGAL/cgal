// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_ROOT_FOR_CIRCLES_2_2_H
#define CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_ROOT_FOR_CIRCLES_2_2_H

#include <CGAL/license/Circular_kernel_2.h>


#include <iostream>
#include <CGAL/Bbox_2.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Handle_for.h>
#include <boost/type_traits/is_same.hpp>

namespace CGAL {

template < typename RT_ >
class Root_for_circles_2_2 {

  typedef RT_                                                              RT;
  typedef typename Root_of_traits< RT >::RootOf_2    Root_of_2;
  typedef typename Root_of_traits< RT >::RootOf_1    FT;

  private:
    Handle_for<Root_of_2> x_;
    Handle_for<Root_of_2> y_;
    
  public:
  Root_for_circles_2_2(){}
  
  Root_for_circles_2_2(const Root_of_2& r1, const Root_of_2& r2)
    : x_(r1), y_(r2)
  {
    // When it is an interval this assertion dont compile
    //CGAL_assertion((r1.is_rational() || r2.is_rational()) || 
    //               (r1.gamma() == r2.gamma()));
  }

  const Root_of_2& x() const 
  { return get_pointee_or_identity(x_); }
    
  const Root_of_2& y() const 
  { return get_pointee_or_identity(y_); }

  CGAL::Bbox_2 bbox() const
  {
    CGAL::Interval_nt<> 
     ix=to_interval(x()),
     iy=to_interval(y());
    return CGAL::Bbox_2(ix.inf(),iy.inf(),
	                ix.sup(),iy.sup());
    /*
      const Root_of_2 &ox = x();
      const Root_of_2 &oy = y();

      if(ox.is_rational() || oy.is_rational()) {
        CGAL::Interval_nt<> 
          ix=to_interval(ox),
          iy=to_interval(oy);
        return CGAL::Bbox_2(ix.inf(),iy.inf(),
	                ix.sup(),iy.sup());
      }

      // delta must be the same
      // WE HAVE TO TEST THE EXECUTION TIME
      // IT STILL NOT POSSIBLE BECAUSE OF THE 
      // PROBLEM ON THE ARRANGEMENT
      // (it is very likely to make it better with this changing)
      const CGAL::Interval_nt<true> alpha1 = to_interval(ox.alpha());
      const CGAL::Interval_nt<true> beta1 = to_interval(ox.beta());
      const CGAL::Interval_nt<true> alpha2 = to_interval(oy.alpha());
      const CGAL::Interval_nt<true> beta2 = to_interval(oy.beta());
      const CGAL::Interval_nt<true> g = to_interval(ox.gamma());
      const CGAL::Interval_nt<true> sqrtg = CGAL::sqrt(g);
      const CGAL::Interval_nt<true> ix = alpha1 + beta1 * sqrtg;
      const CGAL::Interval_nt<true> iy = alpha2 + beta2 * sqrtg;
      return CGAL::Bbox_2(ix.inf(),iy.inf(),
	                ix.sup(),iy.sup());
    */
  }

  template < typename RT >
  friend bool operator == ( const Root_for_circles_2_2<RT>& r1,
   	                    const Root_for_circles_2_2<RT>& r2 );

};
  
template < typename RT >
bool 
operator == ( const Root_for_circles_2_2<RT>& r1,
	      const Root_for_circles_2_2<RT>& r2 )
{ if (CGAL::identical(r1.x_, r2.x_) && CGAL::identical(r1.y_, r2.y_))
	  return true;
  return (r1.x() == r2.x()) && (r1.y() == r2.y()); 
}

template < typename RT >
std::ostream &
operator<<(std::ostream & os, const Root_for_circles_2_2<RT> &r)
{ return os << r.x() << " " << r.y() << " "; }

template < typename RT >
std::istream &
operator>>(std::istream & is, Root_for_circles_2_2<RT> &r)
{
  typedef typename Root_of_traits< RT >::RootOf_2         Root_of_2;
  Root_of_2 x,y;
  
  is >> x >> y;
  if(is)
    r = Root_for_circles_2_2<RT>(x,y);
  return is;
}

} //namespace CGAL

#endif // CGAL_ALGEBRAIC_KERNEL_FOR_CIRCLES_ROOT_FOR_CIRCLES_2_2_H
