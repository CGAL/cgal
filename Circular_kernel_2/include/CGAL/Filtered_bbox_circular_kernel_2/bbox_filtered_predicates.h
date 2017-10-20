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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Constantinos Tsirogiannis , Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)


#ifndef CGAL_BBOX_FILTERED_PREDICATES_H  
#define CGAL_BBOX_FILTERED_PREDICATES_H  

#include <CGAL/license/Circular_kernel_2.h>


#include <CGAL/assertions.h>
#include <CGAL/enum.h>
#include <CGAL/Bbox_2.h>

namespace CGAL {

namespace Bbox_functors {

template <class BK>
class Compare_x_2 : public BK::Circular_kernel:: template Base< BK >::Type::Compare_x_2
{
  typedef typename BK::Circular_arc_point_2                    Circular_arc_point_2;
  typedef typename BK::Point_2                                 Point_2;
  typedef typename BK::Circular_kernel:: 
    template Base< BK >::Type::Compare_x_2                     CK_Compare_x_2;
  typedef CK_Compare_x_2 Base;

public:

  typedef typename CK_Compare_x_2::result_type result_type;

#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else
 template <typename T1, typename T2>
  result_type
  operator()(const T1& t1, const T2& t2) const
  {
    return Base()(t1,t2);
  }
  template <typename T1, typename T2, typename T3>
  result_type
  operator()(const T1& t1, const T2& t2, const T3& t3) const
  {
    return Base()(t1,t2,t3);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  result_type
  operator()(const T1& t1, const T2& t2, const T3& t3, const T4& t4) const
  {
    return Base()(t1,t2,t3,t4);
  }

#endif
  
	
  result_type
  operator()( const Circular_arc_point_2 &a, const Circular_arc_point_2 &b) const
  {
    Bbox_2 bb1=a.bbox(),bb2=b.bbox();

    if( bb1.xmin()>bb2.xmax() )
      return LARGER;

    if( bb1.xmax()<bb2.xmin() )
      return SMALLER;

    return CK_Compare_x_2()(a,b);
  }

};


template <class BK>
class Compare_y_2 : public BK::Circular_kernel:: template Base< BK >::Type::Compare_y_2
{
  typedef typename BK::Circular_arc_point_2                    Circular_arc_point_2;
  typedef typename BK::Point_2                                 Point_2;
  typedef typename BK::Circular_kernel:: 
    template Base< BK >::Type::Compare_y_2                     CK_Compare_y_2;
  typedef CK_Compare_y_2 Base;
public:

  typedef typename CK_Compare_y_2::result_type result_type;

  result_type
  operator() (const Point_2 &p0,
              const Point_2 &p1) const
  {
    return CK_Compare_y_2()(p0, p1);
  }

#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else
  template <typename T1, typename T2>
  result_type
  operator()(const T1& t1, const T2& t2) const
  {
    return Base()(t1,t2);
  }
  template <typename T1, typename T2, typename T3>
  result_type
  operator()(const T1& t1, const T2& t2, const T3& t3) const
  {
    return Base()(t1,t2,t3);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  result_type
  operator()(const T1& t1, const T2& t2, const T3& t3, const T4& t4) const
  {
    return Base()(t1,t2,t3,t4);
  }
#endif
  
  result_type
  operator()( const Circular_arc_point_2 &a, const Circular_arc_point_2 &b) const
  {
    Bbox_2 bb1=a.bbox(),bb2=b.bbox();

    if( bb1.ymin()>bb2.ymax() )
      return LARGER;

    if( bb1.ymax()<bb2.ymin() )
      return SMALLER;

    return CK_Compare_y_2()(a,b);
  }

};

template <class BK>
class Compare_xy_2 : public BK::Circular_kernel:: template Base< BK >::Type::Compare_xy_2
{
  typedef typename BK::Circular_kernel::
  template Base< BK >::Type::Compare_xy_2                      CK_Compare_xy_2;
  typedef CK_Compare_xy_2 Base;
  typedef typename BK::Circular_arc_point_2                    Circular_arc_point_2;
  typedef typename BK::Point_2                                 Point_2;

public:

  typedef typename Base::result_type result_type;


#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else
 template <typename T1, typename T2>
  result_type
  operator()(const T1& t1, const T2& t2) const
  {
    return Base()(t1,t2);
  }
#endif
public:

  result_type
  operator()( const Circular_arc_point_2 &a, const Circular_arc_point_2 &b) const
  {
    typename BK::Compare_x_2 compx;
    typename BK::Compare_y_2 compy;
    Comparison_result tmp;

    if( (tmp=compx(a,b))!=EQUAL)
      return tmp;

    return compy(a,b);
  }

};


template <class BK>
class In_x_range_2 : public BK::Circular_kernel:: template Base< BK >::Type::In_x_range_2
{
  typedef typename BK::Circular_kernel:: 
    template Base< BK >::Type::In_x_range_2                    CK_In_x_range_2;
  typedef CK_In_x_range_2 Base;
  typedef typename BK::Circular_arc_point_2                    Circular_arc_point_2;
  typedef typename BK::Circular_arc_2                          Circular_arc_2;
  typedef typename BK::Line_arc_2                              Line_arc_2;

public:

  typedef typename CK_In_x_range_2::result_type result_type;

#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else
 template <typename T1, typename T2>
  result_type
  operator()(const T1& t1, const T2& t2) const
  {
    return Base()(t1,t2);
  }
#endif
private:

  template <class Arc_2>
  result_type
  _in_x_range_2(const Arc_2 &a, const Circular_arc_point_2 &p) const 
  {

    Bbox_2 bb11 = a.source().bbox(),
      bb12 = a.target().bbox(),
      bb2=p.bbox();      
    if(bb11.xmin() > bb12.xmax()) {
      if(bb2.xmax() < bb12.xmin()) return false;
      else if(bb2.xmin() > bb11.xmax()) return false;
      else if(bb12.xmax() < bb2.xmin() &&
              bb2.xmax() < bb11.xmin()) return true;
    } else if(bb11.xmax() < bb12.xmin()) {
      if(bb2.xmax() < bb11.xmin()) return false;
      else if(bb2.xmin() > bb12.xmax()) return false;
      else if(bb11.xmax() < bb2.xmin() &&
              bb2.xmax() < bb12.xmin()) return true;
    } else {
      if(bb2.xmin() > (std::max)(bb11.xmax(),bb12.xmax())) return false;
      if(bb2.xmax() < (std::min)(bb11.xmin(),bb12.xmin())) return false;
    }
    
    CK_In_x_range_2 Range;

    return Range(a,p);

  }

public:

  result_type
  operator()( const Circular_arc_2 &a, const Circular_arc_point_2 &p) const
  { 
    CGAL_precondition( a.is_x_monotone());
    return _in_x_range_2(a,p);
  }

  result_type
  operator()( const Line_arc_2 &a, const Circular_arc_point_2 &p) const
  { return _in_x_range_2(a,p);}
    

};


template <class BK>
class Compare_y_at_x_2 : public BK::Circular_kernel:: template Base< BK >::Type::Compare_y_at_x_2
{
  typedef typename BK::Circular_kernel:: 
    template Base< BK >::Type::Compare_y_at_x_2                         CK_Compare_y_at_x_2;
  typedef CK_Compare_y_at_x_2 Base;
  typedef typename BK::Circular_arc_2                                   Circular_arc_2;
  typedef typename BK::Circular_arc_point_2                             Circular_arc_point_2;
  typedef typename BK::Line_arc_2                                       Line_arc_2;

public:

  typedef typename CK_Compare_y_at_x_2::result_type result_type;

#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else
  template <typename T1, typename T2>
  result_type
  operator()(const T1& t1, const T2& t2) const
  {
    return Base()(t1,t2);
  } 

  template <typename T1, typename T2, typename T3>
  result_type
  operator()(const T1& t1, const T2& t2, const T3& t3) const
  {
    return Base()(t1,t2,t3);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  result_type
  operator()(const T1& t1, const T2& t2, const T3& t3, const T4& t4) const
  {
    return Base()(t1,t2,t3,t4);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  result_type
  operator()(const T1& t1, const T2& t2, const T3& t3, const T4& t4, const T5& t5) const
  {
    return Base()(t1,t2,t3,t4,t5);
  }
#endif
private:

  template <class Arc_2>
  result_type
  _compare_y_at_x_2(const Circular_arc_point_2 &p,const Arc_2 &a) const
  {
    CGAL_precondition_code(bool tmp=In_x_range_2<BK>()(a,p));
    CGAL_precondition(tmp );

    Bbox_2 bb1=a.bbox(),bb2=p.bbox();

    if(bb1.ymin()>bb2.ymax())
      return SMALLER;
    else if(bb1.ymax()<bb2.ymin())
      return LARGER;

    return CK_Compare_y_at_x_2()(p,a);

  }

public:

  result_type
  operator()( const Circular_arc_point_2 &p,const Circular_arc_2 &a ) const
  {   
    CGAL_precondition( a.is_x_monotone());
    return _compare_y_at_x_2(p,a);
  }

  result_type
  operator()( const Circular_arc_point_2 &p,const Line_arc_2 &a ) const
  {return _compare_y_at_x_2(p,a);}

};


template <class BK>
class Has_on_2 : public BK::Circular_kernel:: template Base< BK >::Type::Has_on_2
{
  typedef typename BK::Circular_kernel:: 
    template Base< BK >::Type::Has_on_2                                 CK_Has_on_2;
  typedef CK_Has_on_2 Base;
  typedef typename BK::Circular_arc_2                                   Circular_arc_2;
  typedef typename BK::Circular_arc_point_2                             Circular_arc_point_2;
  typedef typename BK::Line_arc_2                                       Line_arc_2;

public:

  typedef typename CK_Has_on_2::result_type result_type;
  #ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else
 template <typename T1, typename T2>
  result_type
  operator()(const T1& t1, const T2& t2) const
  {
    return Base()(t1,t2);
  }
#endif
private:

  template <class Arc_2>
  result_type
  _has_on_2(const Arc_2 &a, const Circular_arc_point_2 &p) const
  {
    Bbox_2 bb1=a.bbox(),bb2=p.bbox();

    if(do_overlap(bb1,bb2))
      return CK_Has_on_2()(a,p);

    return false;
  }

public:

  result_type
  operator()( const Circular_arc_2 &a,const Circular_arc_point_2 &p ) const
  {     
    CGAL_precondition( a.is_x_monotone());
    return _has_on_2(a,p);
  }

  result_type
  operator()( const Line_arc_2 &a, const Circular_arc_point_2 &p ) const
  {return _has_on_2(a,p);}

};


template <class BK>
class Equal_2
#ifndef CGAL_CFG_MATCHING_BUG_6
  : public BK::Circular_kernel:: template Base< BK >::Type::Equal_2
#endif
{
  typedef typename BK::Circular_kernel:: 
    template Base< BK >::Type::Equal_2                                  CK_Equal_2;
  typedef typename BK::Circular_arc_2                                   Circular_arc_2;
  typedef typename BK::Point_2                                          Point_2;
  typedef typename BK::Direction_2                                      Direction_2;
  typedef typename BK::Vector_2                                         Vector_2;
  typedef typename BK::Segment_2                                        Segment_2 ;
  typedef typename BK::Ray_2                                            Ray_2;
  typedef typename BK::Line_2                                           Line_2;
  typedef typename BK::Circle_2                                         Circle_2;
  typedef typename BK::Triangle_2                                       Triangle_2;
  typedef typename BK::Iso_rectangle_2                                  Iso_rectangle_2;
  typedef typename BK::Circular_arc_point_2                             Circular_arc_point_2;
  typedef typename BK::Line_arc_2                                       Line_arc_2;

  typedef CK_Equal_2 Base;

public:

  typedef typename CK_Equal_2::result_type result_type;

#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else  
 template <typename T1, typename T2>
  result_type
  operator()(const T1& t1, const T2& t2) const
  {
    return Base()(t1,t2);
  }
#endif

private:

  template <class Arc_2>
  result_type
  _equal_2(const Arc_2 &a,const Arc_2 &b) const
  {
    Bbox_2 bb11=a.source().bbox(),
      bb12=a.target().bbox(),
      bb21=b.source().bbox(),
      bb22=b.target().bbox();

    if(bb11.xmin() > bb21.xmax()) return false;
    if(bb11.xmax() < bb21.xmin()) return false;
    if(bb11.ymin() > bb21.ymax()) return false;
    if(bb11.ymax() < bb21.ymin()) return false;

    if(bb12.xmin() > bb22.xmax()) return false;
    if(bb12.xmax() < bb22.xmin()) return false;
    if(bb12.ymin() > bb22.ymax()) return false;
    if(bb12.ymax() < bb22.ymin()) return false;

    return CK_Equal_2()( a,b );

  }

public:

  result_type
  operator()( const Circular_arc_point_2 &a ,
              const Circular_arc_point_2 &b) const
  { 
    Bbox_2 bb1=a.bbox(),bb2=b.bbox();
    if(bb1.xmin() > bb2.xmax()) return false;
    if(bb1.xmax() < bb2.xmin()) return false;
    if(bb1.ymin() > bb2.ymax()) return false;
    if(bb1.ymax() < bb2.ymin()) return false;
    return CK_Equal_2()( a,b );
  }

  /* WAS THAT HERE FOR OTHER COMPILERS THAN VC* ???
  // redefine to solve ambiguous call error
  result_type
  operator()( const Point_2 &a ,
  const Point_2 &b) const
  { 
  return CK_Equal_2()( a, b);
  }
  */
  result_type
  operator()( const Circular_arc_2 &a , const Circular_arc_2 &b ) const
  {
    CGAL_precondition( a.is_x_monotone());
    CGAL_precondition( b.is_x_monotone());

    return _equal_2(a,b);      
  }

  result_type
  operator()( const Line_arc_2 &a ,
              const Line_arc_2 &b ) const
  {  return _equal_2(a,b);}

  result_type
  operator()( const Circular_arc_2 & ,
              const Line_arc_2 & ) const
  {  return false;}

  result_type
  operator()( const Line_arc_2 & ,
              const Circular_arc_2 & ) const
  {  return false;}

};


template <class BK>
class Do_overlap_2 : public BK::Circular_kernel:: template Base< BK >::Type::Do_overlap_2
{
  typedef typename BK::Circular_kernel:: 
    template Base< BK >::Type::Do_overlap_2                           CK_Do_overlap_2;
  typedef CK_Do_overlap_2 Base;
  typedef typename BK::Circular_arc_2                                 Circular_arc_2;
  typedef typename BK::Line_arc_2                                     Line_arc_2;

public:

  typedef typename CK_Do_overlap_2::result_type result_type;

#ifndef CGAL_CFG_MATCHING_BUG_6
  using Base::operator();
#else
 template <typename T1, typename T2>
  result_type
  operator()(const T1& t1, const T2& t2) const
  {
    return Base()(t1,t2);
  }
#endif

private:

  template <class Arc_2>
  result_type
  _do_overlap_2(const Arc_2 &a, const Arc_2 &b) const
  {
    Bbox_2 bb1=a.bbox(),bb2=b.bbox();  

    if(do_overlap(bb1,bb2))
      return CK_Do_overlap_2()(a,b);
        
    return false;        
  }


public:
    
  result_type
  operator()( const Circular_arc_2 &a , const Circular_arc_2 &b ) const
  {
    CGAL_precondition( a.is_x_monotone());
    CGAL_precondition( b.is_x_monotone());
    return _do_overlap_2(a,b); 
  }

  result_type
  operator()( const Line_arc_2 &a ,
              const Line_arc_2 &b ) const
  {  return _do_overlap_2(a,b);}

  result_type
  operator()( const Circular_arc_2 & ,
              const Line_arc_2 & ) const
  {  return false;}

  result_type
  operator()( const Line_arc_2 & ,
              const Circular_arc_2 & ) const
  {  return false;}

};


template < class BK >
class Intersect_2 : public BK::Circular_kernel:: template Base< BK >::Type::Intersect_2
{
public:
  typedef typename BK::Circular_kernel:: 
    template Base< BK >::Type::Intersect_2      CK_Intersect_2;

  typedef typename BK::Circular_arc_2           Circular_arc_2;
  typedef typename BK::Circular_arc_point_2     Circular_arc_point_2;
  typedef typename BK::Line_arc_2               Line_arc_2;
  typedef typename BK::Circle_2                 Circle;
  typedef typename BK::Line_2                   Line_2;

  using CK_Intersect_2::operator();

  template < class OutputIterator >
  OutputIterator
  operator()(const Line_2 & c1, const Circle & c2, OutputIterator res)
  {
    return CK_Intersect_2()(c1,c2,res);
  }
	
  template < class OutputIterator >
  OutputIterator
  operator()(const Circle & c1, const Line_2 & c2, OutputIterator res)
  {
    return CK_Intersect_2()(c1,c2,res); 
  }

  template < class OutputIterator >
  OutputIterator
  operator()(const Line_arc_2 & c1, const Circle & c2, OutputIterator res)
  {
    if(!do_overlap(c1.bbox(),c2.bbox()))
      return res;

    return CK_Intersect_2()(c1,c2,res);
  }

  template < class OutputIterator >
  OutputIterator
  operator()(const Circle & c1, const Line_arc_2 & c2, OutputIterator res)
  {
    if(!do_overlap(c1.bbox(),c2.bbox()))
      return res;
    return CK_Intersect_2()(c1,c2,res); 
  }

  template < class OutputIterator >
  OutputIterator
  operator()(const Line_2 & c1, const Circular_arc_2 & c2, 
             OutputIterator res)
  {
    return CK_Intersect_2()(c1,c2,res); 
  }

  template < class OutputIterator >
  OutputIterator
  operator()(const Line_2 & c1, const Line_arc_2 & c2, 
             OutputIterator res)
  {
    return CK_Intersect_2()(c1,c2,res); 
  }

  template < class OutputIterator >
  OutputIterator
  operator()(const Circle & c1, const Circle & c2, OutputIterator res)
  { 
    if(!do_overlap(c1.bbox(),c2.bbox()))
      return res;

    return CK_Intersect_2()(c1,c2,res); 
  }

  template < class OutputIterator >
  OutputIterator
  operator()(const Circle & c1, const Circular_arc_2 & c2, OutputIterator res)
  { 
    return operator()(Circular_arc(c1),c2,res);
  }

  template < class OutputIterator >
  OutputIterator
  operator()(const Circular_arc_2 & c1, const Circle & c2, OutputIterator res)
  { 
    return operator()(c1,Circular_arc_2(c2),res);
  }

  template < class OutputIterator >
  OutputIterator
  operator()(const Circular_arc_2 & c1, const Circular_arc_2 & c2, 
             OutputIterator res)
  { 
    Bbox_2 bb1=c1.bbox(),bb2=c2.bbox();

    if(!do_overlap(bb1,bb2 ))
      return res;
	
    return CK_Intersect_2()(c1,c2,res); 
  }

  template < class OutputIterator >
  OutputIterator
  operator()(const Line_arc_2 & c1, const Line_arc_2 & c2, 
             OutputIterator res)
  { 
    Bbox_2 bb1=c1.bbox(),bb2=c2.bbox();

    if(!do_overlap(bb1,bb2 ))
      return res;

    return CK_Intersect_2()(c1,c2,res); 
  }

  template < class OutputIterator >
  OutputIterator
  operator()(const Circular_arc_2 & c1, const Line_arc_2 & c2, 
             OutputIterator res)
  { 
    Bbox_2 bb1=c1.bbox(),bb2=c2.bbox();

    if(!do_overlap(bb1,bb2 ))
      return res;

    return CK_Intersect_2()(c1,c2,res); 
  }

  template < class OutputIterator >
  OutputIterator
  operator()(const Line_arc_2 & c1, const Circular_arc_2 & c2, 
             OutputIterator res)
  {	return operator()(c2,c1,res);}

  template < class OutputIterator >
  OutputIterator
  operator()(const Circular_arc_2 & c1, const Line_2 & c2, 
             OutputIterator res)
  {	return operator()(c2,c1,res);}

  template < class OutputIterator >
  OutputIterator
  operator()(const Line_arc_2 & c1, const Line_2 & c2, 
             OutputIterator res)
  {	return operator()(c2,c1,res);}
   
};

  
} //Bbox_functors

} //namespace CGAL  

#endif // CGAL_BBOX_FILTERED_PREDICATES_H  
