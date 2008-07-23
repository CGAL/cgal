// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Constantinos Tsirogiannis , Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_CIRCULAR_ARC_POINT_WITH_BBOX_H 
#define CGAL_CIRCULAR_ARC_POINT_WITH_BBOX_H 

#include<CGAL/Bbox_2.h>

CGAL_BEGIN_NAMESPACE

template < class BK>
class Circular_arc_point_with_bbox_2 {

    typedef typename BK::Circular_kernel                         CK;
    typedef typename CK::FT                                    FT;
    typedef typename CK::RT                                    RT;
    typedef typename CK::Point_2                               Point_2;
    typedef typename CK::Line_2                                Line_2;
    typedef typename CK::Circle_2                              Circle_2;
    typedef typename CK::Circular_arc_point_2                  RCircular_arc_point_2;
    typedef typename CK::Circular_arc_2                        Circular_arc_2;
    typedef typename CK::Root_of_2                             Root_of_2;

public:
    typedef typename RCircular_arc_point_2::Root_for_circles_2_2 
     Root_for_circles_2_2;
    typedef CK   R; 
  


  ////Construction/////
  Circular_arc_point_with_bbox_2()
    : P_point(),bb(NULL)
    {}

  Circular_arc_point_with_bbox_2(const Root_for_circles_2_2 & np)
    : P_point(np),bb(NULL)
      {}

  Circular_arc_point_with_bbox_2(const RCircular_arc_point_2 & p)
    : P_point(p),bb(NULL)
      {}

  // This avoids Memory Leaks, but may decrease the performance
  // probably not the best solution
  Circular_arc_point_with_bbox_2(const Circular_arc_point_with_bbox_2 &c) : P_point(c.P_point), bb(NULL) { }
	~Circular_arc_point_with_bbox_2() { if(bb) delete bb; }

  ////Accesors////
  const RCircular_arc_point_2 & point() const
  {return P_point;}
            
  typename Qualified_result_of<typename R::Compute_Circular_x_2,RCircular_arc_point_2>::type
  x() const
    { return P_point.x();}

  typename Qualified_result_of<typename R::Compute_Circular_y_2,RCircular_arc_point_2>::type
  y() const
    { return P_point.y();}


  ////Bbox related accessors////
  
bool has_no_bbox() const
  { return (bb==NULL);}

  Bbox_2  bbox() const
    { 
      if(this->has_no_bbox())
        bb= new Bbox_2(P_point.bbox());
              
        return *bb;     
    }


private:

   RCircular_arc_point_2  P_point;
   mutable Bbox_2         *bb;

};

template < typename BboxKernel >
inline
bool
operator==(const Circular_arc_point_with_bbox_2<BboxKernel> &p,
           const Circular_arc_point_with_bbox_2<BboxKernel> &q)
{
  return BboxKernel().equal_2_object()(p, q);
}

template < typename  BboxKernel >
inline
bool
operator!=(const Circular_arc_point_with_bbox_2<BboxKernel> &p,
           const Circular_arc_point_with_bbox_2<BboxKernel> &q)
{
  return ! (p == q);
}

  template < typename BK >
  std::istream &
  operator>>(std::istream & is, Circular_arc_point_with_bbox_2<BK> &p)
  {
    typedef typename BK::CK                      CK;
    typedef typename CK::Root_of_2               Root_of_2;
    typedef typename CK::Root_for_circles_2_2    Root_for_circles_2_2;

    Root_for_circles_2_2 r;
    is >> r;
    if(is)
      p = Circular_arc_point_with_bbox_2<BK>(r);
    return is;
  }

template < class BK >
std::ostream&
operator<<(std::ostream &os, const Circular_arc_point_with_bbox_2<BK> &p)
{
  //I can make it because I know the output format of Root_for_circle
    return os << p.x() << " " << p.y() << " ";
}

CGAL_END_NAMESPACE

#endif // CGAL_CIRCULAR_ARC_POINT_WITH_BBOX_H 
