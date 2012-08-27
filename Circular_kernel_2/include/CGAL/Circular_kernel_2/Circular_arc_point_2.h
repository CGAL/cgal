// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_CIRCULAR_ARC_POINT_2_H
#define CGAL_CIRCULAR_KERNEL_CIRCULAR_ARC_POINT_2_H

#include <iostream>
#include <CGAL/Bbox_2.h>
#include <CGAL/Interval_nt.h>
#include <boost/type_traits/is_same.hpp>

namespace CGAL {
namespace internal {

  template <class CK >
  class Circular_arc_point_2_base
  {
    typedef typename CK::FT                      FT;
    typedef typename CK::Root_of_2               Root_of_2;
    typedef typename CK::Point_2                 Point_2;
    
  public: // fixme ?
    typedef typename CK::Root_for_circles_2_2 Root_for_circles_2_2;
    typedef typename CK::template Handle<Root_for_circles_2_2>::type  Base;
    
    Circular_arc_point_2_base() 
    {}
    
    Circular_arc_point_2_base(const Root_for_circles_2_2 & np)
      :  _p(np)
    {}

    Circular_arc_point_2_base(const Point_2 & p)
      :  _p(p.x(),p.y()/*,1,1,-p.x()-p.y()*/)
    {}

    const Root_of_2 & x() const 
    { return get(_p).x(); }
    
    const Root_of_2 & y() const 
    { return get(_p).y(); }
    
    CGAL::Bbox_2 bbox() const
    {
      return get(_p).bbox();
    }

    const Root_for_circles_2_2 & coordinates() const 
    { return get(_p); }

    bool equal_ref(const Circular_arc_point_2_base &p) const
    {
      return CGAL::identical(_p, p._p);      
    }

  private:
    Base _p;
  };

  template < typename CK >
  std::ostream &
  print(std::ostream & os, const Circular_arc_point_2_base<CK> &p)
  {
    return os << "CirclArcEndPoint_2(" << std::endl
	      << p.x() << ", " << p.y() << ')';
  }

template < typename BK, typename Base_CK >
class Filtered_bbox_circular_arc_point_2_base
  : public Base_CK::Circular_arc_point_2
{
public:
  typedef Filtered_bbox_circular_arc_point_2_base<BK,Base_CK> Self;
  typedef typename Base_CK::Circular_arc_point_2 P_point;

  typedef typename BK::Point_2               Point_2;
  typedef typename BK::Root_for_circles_2_2  Root_for_circles_2_2;

  ////Construction/////
  Filtered_bbox_circular_arc_point_2_base()
    : P_point(), bb(NULL)
  {}

  Filtered_bbox_circular_arc_point_2_base(const P_point& pt)
    : P_point(pt), bb(NULL)
  {}

  explicit Filtered_bbox_circular_arc_point_2_base(const Root_for_circles_2_2 & np)
    : P_point(np), bb(NULL)
  {}

  explicit Filtered_bbox_circular_arc_point_2_base(const Point_2 & p)
    : P_point(p), bb(NULL)
  {}

  Filtered_bbox_circular_arc_point_2_base(const Self &c) 
    : P_point(c), bb(c.bb ? new Bbox_2(*(c.bb)) : NULL)
  {}

  Filtered_bbox_circular_arc_point_2_base&
  operator=(const Self& c) {
    if(this != &c)
    {
      this->P_point::operator=(c);
      bb = c.bb ? new Bbox_2(*(c.bb)) : NULL;
    }
    return *this;
  }
	
  ~Filtered_bbox_circular_arc_point_2_base() {
    if(bb) {
      delete bb; 
      bb = 0;
    }
  }

  ////Bbox related accessors////
  
  bool has_no_bbox() const
  { return (bb==NULL);}

  Bbox_2  bbox() const
  { 
    if(this->has_no_bbox())
      bb= new Bbox_2(P_point::bbox());
              
    return *bb;     
  }

private:
  mutable Bbox_2         *bb;
}; // end class Filtered_bbox_circular_arc_point_2_base


} // namespace internal
} // namespace CGAL

#endif // CGAL_CIRCULAR_KERNEL_CIRCULAR_ARC_POINT_2_H
