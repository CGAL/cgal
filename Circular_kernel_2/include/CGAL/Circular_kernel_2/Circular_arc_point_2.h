// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_KERNEL_CIRCULAR_ARC_POINT_2_H
#define CGAL_CIRCULAR_KERNEL_CIRCULAR_ARC_POINT_2_H

#include <CGAL/license/Circular_kernel_2.h>


#include <iostream>
#include <CGAL/Handle.h>
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
    { return get_pointee_or_identity(_p).x(); }
    
    const Root_of_2 & y() const 
    { return get_pointee_or_identity(_p).y(); }
    
    CGAL::Bbox_2 bbox() const
    {
      return get_pointee_or_identity(_p).bbox();
    }

    const Root_for_circles_2_2 & coordinates() const 
    { return get_pointee_or_identity(_p); }

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
    : P_point(), bb(nullptr)
  {}

  Filtered_bbox_circular_arc_point_2_base(const P_point& pt)
    : P_point(pt), bb(nullptr)
  {}

  explicit Filtered_bbox_circular_arc_point_2_base(const Root_for_circles_2_2 & np)
    : P_point(np), bb(nullptr)
  {}

  explicit Filtered_bbox_circular_arc_point_2_base(const Point_2 & p)
    : P_point(p), bb(nullptr)
  {}

  Filtered_bbox_circular_arc_point_2_base(const Self &c) 
    : P_point(c), bb(c.bb ? new Bbox_2(*(c.bb)) : nullptr)
  {}

  Filtered_bbox_circular_arc_point_2_base&
  operator=(const Self& c) {
    if(this != &c)
    {
      this->P_point::operator=(c);

      if (bb != nullptr){ 
        delete bb;
      }
      bb = c.bb ? new Bbox_2(*(c.bb)) : nullptr;
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
  { return (bb==nullptr);}

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
