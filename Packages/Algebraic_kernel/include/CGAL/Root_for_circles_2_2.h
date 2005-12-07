// Copyright (c) 2005  INRIA Sophia-Antipolis (France) 
// All rights reserved.
//
// Authors : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//           Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
// 
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_ROOT_FOR_CIRCLES_2_2_H
#define CGAL_ROOT_FOR_CIRCLES_2_2_H

#include <iostream>

CGAL_BEGIN_NAMESPACE

template < typename RT_ >
class Root_for_circles_2_2 {

  typedef RT_                                                              RT;
  typedef typename Root_of_traits< RT >::RootOf_2    Root_of_2;

  private:
    Root_of_2 x_;
    Root_of_2 y_;
    
  public:
  Root_for_circles_2_2(){}
    
  Root_for_circles_2_2(const Root_of_2& r1, const Root_of_2& r2)
    : x_(r1), y_(r2)
  {}

  const Root_of_2& x() const 
  { return x_; }
    
  const Root_of_2& y() const 
  { return y_; }
};
  
template < typename RT >
bool 
operator == ( const Root_for_circles_2_2<RT>& r1,
	      const Root_for_circles_2_2<RT>& r2 )
{ return (r1.x() == r2.x()) && (r1.y() == r2.y()); }

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
