// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud, Sylvain Pion

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ROOT_FOR_CIRCLES_2_2_H
#define CGAL_ROOT_FOR_CIRCLES_2_2_H

#include <iostream>
#include <CGAL/Bbox_2.h>
#include <CGAL/Handle_for.h>

CGAL_BEGIN_NAMESPACE

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
  {}

  const Root_of_2& x() const 
  { return get(x_); }
    
  const Root_of_2& y() const 
  { return get(y_); }

  CGAL::Bbox_2 bbox() const
  {
    CGAL::Interval_nt<> 
        ix=to_interval(x()),
        iy=to_interval(y());
    return CGAL::Bbox_2(ix.inf(),iy.inf(),
	                ix.sup(),iy.sup());
  }

};
  
template < typename RT >
bool 
operator == ( const Root_for_circles_2_2<RT>& r1,
	      const Root_for_circles_2_2<RT>& r2 )
{ if (CGAL::identical(r1.x_, r2.x_) && CGAL::identical(r1.x_, r2.x_))
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

CGAL_END_NAMESPACE

#endif // CGAL_ROOT_FOR_CIRCLES_2_2_H
