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


#ifndef CGAL_BBOX_FILTERED_PREDICATES_H  
#define CGAL_BBOX_FILTERED_PREDICATES_H  

#include <CGAL/assertions.h>
#include <CGAL/enum.h>
#include <CGAL/Object.h>
#include <CGAL/Bbox_2.h>

CGAL_BEGIN_NAMESPACE

namespace Bbox_functors {

template <class BK>
class Compare_x_2 : public BK::Circular_kernel::Compare_x_2
  {
    typedef typename BK::Circular_kernel                         CK;
    typedef typename BK::Circular_arc_point_2                    Circular_arc_point_2;

   public:

    typedef typename CK::Compare_x_2::result_type result_type; 
    using CK::Compare_x_2::operator();

   public:
	
    result_type
    operator()( const Circular_arc_point_2 &a, const Circular_arc_point_2 &b) const
    {
       Bbox_2 bb1=a.bbox(),bb2=b.bbox();

       if( bb1.xmin()>bb2.xmax() )
         return LARGER;

       if( bb1.xmax()<bb2.xmin() )
         return SMALLER;

       return CK().compare_x_2_object()(a.point(),b.point());
    }

  };


template <class BK>
class Compare_y_2 : public BK::Circular_kernel::Compare_y_2
  {
    typedef typename BK::Circular_kernel                         CK;
    typedef typename BK::Circular_arc_point_2                    Circular_arc_point_2;

   public:

    typedef typename CK::Compare_y_2::result_type result_type; 
    using CK::Compare_y_2::operator();

   public:

    result_type
    operator()( const Circular_arc_point_2 &a, const Circular_arc_point_2 &b) const
    {
       Bbox_2 bb1=a.bbox(),bb2=b.bbox();


       if( bb1.ymin()>bb2.ymax() )
         return LARGER;

       if( bb1.ymax()<bb2.ymin() )
         return SMALLER;

       return CK().compare_y_2_object()(a.point(),b.point());
    }

  };

template <class BK>
class Compare_xy_2 : public BK::Circular_kernel::Compare_xy_2
  {
    typedef typename BK::Circular_kernel                         CK;
    typedef typename BK::Circular_arc_point_2                    Circular_arc_point_2;

   public:

    typedef typename CK::Compare_xy_2::result_type result_type; 
    using CK::Compare_xy_2::operator();

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
class In_x_range_2 : public BK::Circular_kernel::In_x_range_2
  {
    typedef typename BK::Circular_kernel                         CK;
    typedef typename BK::Circular_arc_point_2                    Circular_arc_point_2;
    typedef typename BK::Circular_arc_2                          Circular_arc_2;
    typedef typename BK::Line_arc_2                              Line_arc_2;

   public:

    typedef typename CK::In_x_range_2::result_type result_type; 
    using CK::In_x_range_2::operator();

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
    
      typename CK::In_x_range_2 Range;

      return Range(a.arc(),p.point());

    }

   public:

    result_type
    operator()( const Circular_arc_2 &a, const Circular_arc_point_2 &p) const
    { 
      CGAL_precondition( a.arc().is_x_monotone());
      return _in_x_range_2(a,p);
   }

    result_type
    operator()( const Line_arc_2 &a, const Circular_arc_point_2 &p) const
    { return _in_x_range_2(a,p);}
    

  };



template < class BK >
  class Construct_line_arc_2
  {

    typedef typename BK::Point_2                   Point_2;
    typedef typename BK::Line_2                    Line_2;
    typedef typename BK::Circle_2                  Circle_2;
    typedef typename BK::Circular_arc_point_2      Circular_arc_point_2;
    typedef typename BK::Segment_2                 Segment_2;
    typedef typename BK::Line_arc_2                Line_arc_2;

  public:
    typedef Line_arc_2   result_type;
    
    result_type
    operator()(void) 
    { return Line_arc_2(); }

    result_type
    operator()(const Line_2 &support,
	       const Circle_2 &c1,const bool b1,
	       const Circle_2 &c2,const bool b2) const
    { return Line_arc_2(support,c1,b1,c2,b2); }

    result_type
    operator()(const Line_2 &support,
	       const Line_2 &l1,
	       const Line_2 &l2) const
    { return Line_arc_2(support,l1,l2); }

    result_type
    operator()(const Line_2 &support,
	       const Circular_arc_point_2 &p1,
	       const Circular_arc_point_2 &p2) const
    { return Line_arc_2(support,p1,p2); }

    result_type
    operator()(const Segment_2 &s) const
    { return Line_arc_2(s); }

    result_type
    operator()(const Point_2 &p1,
	       const Point_2 &p2) const
    { return Line_arc_2(p1,p2); }

  };



  template < class BK >
  class Construct_circular_arc_2
  {

    typedef typename BK::FT                           FT;
    typedef typename BK::RT                           RT;
    typedef typename BK::Point_2                      Point_2;
    typedef typename BK::Line_2                       Line_2;
    typedef typename BK::Circle_2                     Circle_2;
    typedef typename BK::Circular_arc_2               Circular_arc_2;
    typedef typename BK::Circular_arc_point_2         Circular_arc_point_2;

  public:
    typedef  Circular_arc_2 result_type;
    
    result_type
    operator()(void) 
    { return Circular_arc_2(); }

    result_type
    operator()(const Circle_2 &c) const
    { return Circular_arc_2(c); }

    result_type
    operator()(const Circle_2 &support,
               const Circular_arc_point_2 &source, 
               const Circular_arc_point_2 &target) const
    { return Circular_arc_2(support,source,target); }

    result_type
    operator()(const Circle_2 &support,
               const Line_2 &l1, bool b1,
               const Line_2 &l2, bool b2) const
    { return Circular_arc_2(support,l1,b1,l2,b2); }

    result_type
    operator()(const Circle_2 &c,
               const Circle_2 &c1, bool b_1,
               const Circle_2 &c2, bool b_2) const
    { return Circular_arc_2(c,c1,b_1,c2,b_2); }

    result_type
    operator()(const Circular_arc_2 &A,
               bool b,
               const Circle_2 &ccut, bool b_cut) const
    { return Circular_arc_2(A,b,ccut,b_cut); }

    result_type
    operator()(const Point_2 &begin,
               const Point_2 &middle, 
               const Point_2 &end) const
    { return Circular_arc_2(begin,middle,end); }

    result_type
    operator()(const Point_2 &begin,
               const Point_2 &end,
	       const FT& bulge) const
    { return Circular_arc_2(begin,end,bulge); }

  };

  template < class BK >
  class Construct_circular_arc_point_2
  {
    typedef typename BK::Point_2               Point_2;
    typedef typename BK::Circular_arc_point_2  Circular_arc_point_2;
    typedef typename Circular_arc_point_2::Root_for_circles_2_2  
                                               Root_for_circles_2_2;

  public:
    typedef Circular_arc_point_2 result_type;

    result_type
    operator()(void) 
    { return Circular_arc_point_2(); }

    result_type
    operator()(const Root_for_circles_2_2 & np) const
    { return Circular_arc_point_2(np); }

    result_type
    operator()(const Point_2 & p) const
    { return Circular_arc_point_2(p); }

  };

template <class BK>
class Construct_circular_source_vertex_2
  {
    typedef typename BK::Circular_kernel                       CK;
    typedef typename BK::Circular_arc_point_2                Circular_arc_point_2;
    typedef typename BK::Circular_arc_2                      Circular_arc_2;
    typedef typename BK::Line_arc_2                              Line_arc_2;

   public:

    typedef Circular_arc_point_2    result_type;
    typedef const result_type &     qualified_result_type;

    result_type
    operator()(const Circular_arc_2& a) const
    { return CK().construct_circular_source_vertex_2_object()(a.arc()); }

    result_type
    operator()(const Line_arc_2& a) const
    { return CK().construct_circular_source_vertex_2_object()(a.arc()); }

  };


template <class BK>
class Construct_circular_target_vertex_2
  {
    typedef typename BK::Circular_kernel                       CK;
    typedef typename BK::Circular_arc_point_2                Circular_arc_point_2;
    typedef typename BK::Circular_arc_2                      Circular_arc_2;
    typedef typename BK::Line_arc_2                              Line_arc_2;

   public:

    typedef Circular_arc_point_2 result_type;
    typedef const result_type &     qualified_result_type;

    result_type
    operator()(const Circular_arc_2& a) const
    { return CK().construct_circular_target_vertex_2_object()(a.arc()); }

    result_type
    operator()(const Line_arc_2& a) const
    { return CK().construct_circular_target_vertex_2_object()(a.arc()); }

  };



template <class BK>
class Construct_circular_min_vertex_2
  {
    typedef typename BK::Circular_kernel                       CK;
    typedef typename BK::Circular_arc_point_2                Circular_arc_point_2;
    typedef typename BK::Circular_arc_2                      Circular_arc_2;
    typedef typename BK::Line_arc_2                              Line_arc_2;

   public:

    typedef Circular_arc_point_2 result_type;
    
    result_type
    operator()(const Circular_arc_2& a) const
    { return CK().construct_circular_min_vertex_2_object()(a.arc()); }

    result_type
    operator()(const Line_arc_2& a) const
    { return CK().construct_circular_min_vertex_2_object()(a.arc()); }
    
  };

template <class BK>
class Construct_circular_max_vertex_2
  {
    typedef typename BK::Circular_kernel                           CK;
    typedef typename BK::Circular_arc_point_2                    Circular_arc_point_2;
    typedef typename BK::Circular_arc_2                          Circular_arc_2;
    typedef typename BK::Line_arc_2                              Line_arc_2;

   public:

    typedef Circular_arc_point_2 result_type;

    result_type
    operator()(const Circular_arc_2& a) const
    { return CK().construct_circular_max_vertex_2_object()(a.arc()); }

    result_type
    operator()(const Line_arc_2& a) const
    { return CK().construct_circular_max_vertex_2_object()(a.arc()); }

  };

template <class BK>
class Is_vertical_2 : public BK::Circular_kernel::Is_vertical_2
  {
    typedef typename BK::Circular_kernel                         CK;
    typedef typename BK::Circular_arc_point_2                    Circular_arc_point_2;
    typedef typename BK::Circular_arc_2                          Circular_arc_2;
    typedef typename BK::Line_arc_2                              Line_arc_2;

   public:

    typedef typename CK::Is_vertical_2::result_type result_type; 
    using CK::Is_vertical_2::operator();
    
    result_type
      operator()(const Circular_arc_2& a) const
    { return CK().is_vertical_2_object()(a.arc()); }

    result_type
      operator()(const Line_arc_2& a) const
    { return CK().is_vertical_2_object()(a.arc()); }

  };

template <class BK>
class Compare_y_at_x_2 : public BK::Circular_kernel::Compare_y_at_x_2
  {
    typedef typename BK::Circular_kernel                                  CK;
    typedef typename BK::Circular_arc_2                                   Circular_arc_2;
    typedef typename BK::Circular_arc_point_2                             Circular_arc_point_2;
    typedef typename BK::Line_arc_2                                       Line_arc_2;

  public:

    typedef typename CK::Compare_y_at_x_2::result_type result_type; 
    using CK::Compare_y_at_x_2::operator();

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

      return CK().compare_y_at_x_2_object()(p.point(),a.arc());

    }

  public:

    result_type
    operator()( const Circular_arc_point_2 &p,const Circular_arc_2 &a ) const
    {   
      CGAL_precondition( a.arc().is_x_monotone());
      return _compare_y_at_x_2(p,a);
    }

    result_type
    operator()( const Circular_arc_point_2 &p,const Line_arc_2 &a ) const
     {return _compare_y_at_x_2(p,a);}

};


template <class BK>
class Has_on_2 : public BK::Circular_kernel::Has_on_2
  {
    typedef typename BK::Circular_kernel                                  CK;
    typedef typename BK::Circular_arc_2                                   Circular_arc_2;
    typedef typename BK::Circular_arc_point_2                             Circular_arc_point_2;
    typedef typename BK::Line_arc_2                                       Line_arc_2;

  public:

    typedef typename CK::Has_on_2::result_type result_type; 
    using CK::Has_on_2::operator();

  private:

    template <class Arc_2>
    result_type
    _has_on_2(const Arc_2 &a, const Circular_arc_point_2 &p) const
    {
      Bbox_2 bb1=a.bbox(),bb2=p.bbox();

      if(do_overlap(bb1,bb2))
        return CK().has_on_2_object()(a.arc(),p.point());

       return false;
    }

  public:

    result_type
    operator()( const Circular_arc_2 &a,const Circular_arc_point_2 &p ) const
    {     
      CGAL_precondition( a.arc().is_x_monotone());
      return _has_on_2(a,p);
    }

    result_type
    operator()( const Line_arc_2 &a, const Circular_arc_point_2 &p ) const
     {return _has_on_2(a,p);}

};


template <class BK>
class Equal_2
  {
    typedef typename BK::Circular_kernel                                  CK;
    typedef typename BK::Circular_arc_2                                   Circular_arc_2;
    typedef typename BK::Circular_arc_point_2                             Circular_arc_point_2;
    typedef typename BK::Line_arc_2                                       Line_arc_2;

  public:

    typedef bool result_type; 

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

      return CK().equal_2_object()( a.arc(),b.arc() );

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
      return CK().equal_2_object()( a.point(),b.point() );
    }

    result_type
    operator()( const Circular_arc_2 &a , const Circular_arc_2 &b ) const
    {
      CGAL_precondition( a.arc().is_x_monotone());
      CGAL_precondition( b.arc().is_x_monotone());

      return _equal_2(a,b);      
    }

    result_type
    operator()( const Line_arc_2 &a ,
                const Line_arc_2 &b ) const
    {  return _equal_2(a,b);}

    result_type
    operator()( const Circular_arc_2 &a ,
                const Line_arc_2 &b ) const
	 {  return false;}

    result_type
    operator()( const Line_arc_2 &a ,
                const Circular_arc_2 &b ) const
	 {  return false;}

};


template <class BK>
class Do_overlap_2 : public BK::Circular_kernel::Do_overlap_2
  {
    typedef typename BK::Circular_kernel                                CK;
    typedef typename BK::Circular_arc_2                                 Circular_arc_2;
    typedef typename BK::Line_arc_2                                     Line_arc_2;

  public:

    typedef typename CK::Do_overlap_2::result_type result_type; 
    using CK::Do_overlap_2::operator();

  private:

    template <class Arc_2>
    result_type
    _do_overlap_2(const Arc_2 &a, const Arc_2 &b) const
    {
      Bbox_2 bb1=a.bbox(),bb2=b.bbox();  

      if(do_overlap(bb1,bb2))
        return CK().do_overlap_2_object()(a.arc(),b.arc());
        
      return false;        
    }


  public:
    
    result_type
    operator()( const Circular_arc_2 &a , const Circular_arc_2 &b ) const
    {
      CGAL_precondition( a.arc().is_x_monotone());
      CGAL_precondition( b.arc().is_x_monotone());
      return _do_overlap_2(a,b); 
    }

    result_type
    operator()( const Line_arc_2 &a ,
                const Line_arc_2 &b ) const
    {  return _do_overlap_2(a,b);}

    result_type
    operator()( const Circular_arc_2 &a ,
                const Line_arc_2 &b ) const
	 {  return false;}

    result_type
    operator()( const Line_arc_2 &a ,
                const Circular_arc_2 &b ) const
	 {  return false;}

};


// This predicate cannot be filtered

 template < class BK >
  class Compare_y_to_right_2 : public BK::Circular_kernel::Compare_y_to_right_2
  {
    typedef typename BK::Circular_kernel       CK;
    typedef typename BK::Circular_arc_2        Circular_arc_2;
    typedef typename BK::Circular_arc_point_2  Circular_arc_point_2;

  public:
    typedef typename CK::Compare_y_to_right_2::result_type result_type; 
    using CK::Compare_y_to_right_2::operator();
    
    template <typename T1, typename T2>
    result_type
    operator()(const T1 &a1,
               const T2 &a2,
               const Circular_arc_point_2 &p) const
    { return CK().compare_y_to_right_2_object()(a1.arc(), a2.arc(), p.point()); }
    
  };


  template < class BK >
  class Make_x_monotone_2 : public BK::Circular_kernel::Make_x_monotone_2
  {
    typedef typename BK::Circular_kernel          CK;
    typedef typename CK::Circular_arc_2           Rcirc_arc_2;
    typedef typename BK::Circular_arc_2       	  Circular_arc_2;
    typedef typename BK::Line_arc_2               Line_arc_2;

  public:
    typedef typename CK::Make_x_monotone_2::result_type result_type; 
    using CK::Make_x_monotone_2::operator();

    template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_2 &A, OutputIterator res)
      { 
	std::vector<Object> vec;
	
	CK().make_x_monotone_2_object()(A.arc(), std::back_inserter(vec) );

	for(unsigned i=0; i<vec.size() ; i++)
	{
	  const Rcirc_arc_2 *tmp_arc;
	  tmp_arc=object_cast<Rcirc_arc_2>(&vec.at(i));
	  CGAL_assertion(tmp_arc!=NULL);
	  *res++ = make_object( Circular_arc_2(*tmp_arc) );
	}

	return res;
      }
    
    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_2 &A, OutputIterator res)
      { 
	*res++ = make_object(A);
	return res;
      }
  };

  template < class BK >
  class Do_intersect_2 : public BK::Circular_kernel::Linear_kernel::Do_intersect_2
  {
  public:
    typedef typename BK::Circular_kernel            CK;
    typedef typename BK::Circular_arc_2           Circular_arc_2;
    typedef typename BK::Circular_arc_point_2     Circular_arc_point_2;
    typedef typename BK::Line_arc_2               Line_arc_2;
    typedef typename CK::Circular_arc_2           Rcirc_arc_2;
    typedef typename CK::Line_arc_2               Rline_arc_2;
    typedef typename CK::Circular_arc_point_2     Rcirc_arc_point_2;
    typedef typename BK::Circle_2                 Circle_2;
    typedef typename BK::Line_2                   Line_2;

    typedef typename CK::Do_intersect_2::result_type result_type; 
    using typename CK::Linear_kernel::Do_intersect_2::operator();

    result_type
    operator()(const Circular_arc_2 & c1, const Circular_arc_2 & c2) {
			return CK().do_intersect_2_object()(c1.arc(), c2.arc());
	  }

    result_type
    operator()(const Line_arc_2 & c1, const Line_arc_2 & c2) {
			return CK().do_intersect_2_object()(c1.arc(), c2.arc());
	  }

    result_type
    operator()(const Line_arc_2 & c1, const Circle_2 & c2) {
			return CK().do_intersect_2_object()(c1.arc(), c2);
	  }

    result_type
    operator()(const Circle_2 & c1, const Line_arc_2 & c2) {
			return CK().do_intersect_2_object()(c1, c2.arc());
	  }

    result_type
    operator()(const Circular_arc_2 & c1, const Circle_2 & c2) {
			return CK().do_intersect_2_object()(c1.arc(), c2);
	  }

    result_type
    operator()(const Circle_2 & c1, const Circular_arc_2 & c2) {
			return CK().do_intersect_2_object()(c1, c2.arc());
	  }

    result_type
    operator()(const Line_arc_2 & c1, const Circular_arc_2 & c2) {
			return CK().do_intersect_2_object()(c1.arc(), c2.arc());
	  }

    result_type
    operator()(const Circular_arc_2 & c1, const Line_arc_2 & c2) {
			return CK().do_intersect_2_object()(c1.arc(), c2.arc());
	  }

    result_type
    operator()(const Line_2 & c1, const Line_arc_2 & c2) {
			return CK().do_intersect_2_object()(c1, c2.arc());
	  }

    result_type
    operator()(const Line_arc_2 & c1, const Line_2 & c2) {
			return CK().do_intersect_2_object()(c1.arc(), c2);
	  }

    result_type
    operator()(const Line_2 & c1, const Circular_arc_2 & c2) {
			return CK().do_intersect_2_object()(c1, c2.arc());
	  }

    result_type
    operator()(const Circular_arc_2 & c1, const Line_2 & c2) {
			return CK().do_intersect_2_object()(c1.arc(), c2);
	  }

    result_type
    operator()(const Line_2 & c1, const Circle_2 & c2) {
			return CK().do_intersect_2_object()(c1, c2);
	  }

    result_type
    operator()(const Circle_2 & c1, const Line_2 & c2) {
			return CK().do_intersect_2_object()(c1, c2);
	  }

  };

  template < class BK >
  class Intersect_2 : public BK::Circular_kernel::Linear_kernel::Intersect_2
  {
    public:

    typedef typename BK::Circular_kernel            CK;
    typedef typename BK::Circular_arc_2           Circular_arc_2;
    typedef typename BK::Circular_arc_point_2     Circular_arc_point_2;
    typedef typename BK::Line_arc_2               Line_arc_2;
    typedef typename CK::Circular_arc_2           Rcirc_arc_2;
    typedef typename CK::Line_arc_2               Rline_arc_2;
    typedef typename CK::Circular_arc_point_2     Rcirc_arc_point_2;
    typedef typename BK::Circle_2                 Circle;
    typedef typename BK::Line_2                   Line_2;

    typedef typename CK::Intersect_2::result_type result_type; 
    using typename CK::Linear_kernel::Intersect_2::operator();

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_2 & c1, const Circle & c2, OutputIterator res)
      {
	      std::vector<Object> vec;
         
        CK().intersect_2_object()(c1,c2,std::back_inserter(vec));

        for(unsigned i=0; i<vec.size() ; i++)
        {
          const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;

          if ( (tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(&vec.at(i)))!=NULL )
            *res++ = make_object( std::make_pair(Circular_arc_point_2(tmp_point->first),tmp_point->second));
          else
            *res++=vec.at(i);
        }

        return res;
	    }
	
    template < class OutputIterator >
    OutputIterator
    operator()(const Circle & c1, const Line_2 & c2, OutputIterator res)
      {
	      std::vector<Object> vec;
         
        CK().intersect_2_object()(c1,c2,std::back_inserter(vec));

        for(unsigned i=0; i<vec.size() ; i++)
        {
          const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;

          if ( (tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(&vec.at(i)))!=NULL )
            *res++ = make_object( std::make_pair(Circular_arc_point_2(tmp_point->first),tmp_point->second));
          else
            *res++=vec.at(i);
        }

        return res;
	    }

    template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_2 & c1, const Circle & c2, OutputIterator res)
      {
        if(!do_overlap(c1.bbox(),c2.bbox()))
          return res;

	      std::vector<Object> vec;
         
        CK().intersect_2_object()(c1.arc(),c2,std::back_inserter(vec));

        for(unsigned i=0; i<vec.size() ; i++)
        {
          const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;

          if ( (tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(&vec.at(i)))!=NULL )
            *res++ = make_object( std::make_pair(Circular_arc_point_2(tmp_point->first),tmp_point->second));
          else
            *res++=vec.at(i);
        }

        return res;
	    }

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle & c1, const Line_arc_2 & c2, OutputIterator res)
      {
        if(!do_overlap(c1.bbox(),c2.bbox()))
          return res;

	      std::vector<Object> vec;
         
        CK().intersect_2_object()(c1,c2.arc(),std::back_inserter(vec));

        for(unsigned i=0; i<vec.size() ; i++)
        {
          const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;

          if ( (tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(&vec.at(i)))!=NULL )
            *res++ = make_object( std::make_pair(Circular_arc_point_2(tmp_point->first),tmp_point->second));
          else
            *res++=vec.at(i);
        }

        return res;
	    }

		template < class OutputIterator >
		OutputIterator
		operator()(const Line_2 & c1, const Circular_arc_2 & c2, 
		      OutputIterator res)
		  {
	      std::vector<Object> vec;
         
        CK().intersect_2_object()(c1,c2.arc(),std::back_inserter(vec));

        for(unsigned i=0; i<vec.size() ; i++)
        {
          const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;

          if ( (tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(&vec.at(i)))!=NULL )
            *res++ = make_object( std::make_pair(Circular_arc_point_2(tmp_point->first),tmp_point->second));
          else
            *res++=vec.at(i);
        }

        return res;
		  }

		template < class OutputIterator >
		OutputIterator
		operator()(const Line_2 & c1, const Line_arc_2 & c2, 
		       OutputIterator res)
		  {
	      std::vector<Object> vec;
         
        CK().intersect_2_object()(c1,c2,std::back_inserter(vec));

        for(unsigned i=0; i<vec.size() ; i++)
        {
          const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;

          if ( (tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(&vec.at(i)))!=NULL )
            *res++ = make_object( std::make_pair(Circular_arc_point_2(tmp_point->first),tmp_point->second));
          else
            *res++=vec.at(i);
        }

        return res;
			}

    template < class OutputIterator >
    OutputIterator
    operator()(const Circle & c1, const Circle & c2, OutputIterator res)
      { 
        if(!do_overlap(c1.bbox(),c2.bbox()))
          return res;

	      std::vector<Object> vec;
         
        CK().intersect_2_object()(c1,c2,std::back_inserter(vec));

        for(unsigned i=0; i<vec.size() ; i++)
        {
          const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;

          if ( (tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(&vec.at(i)))!=NULL )
            *res++ = make_object( std::make_pair(Circular_arc_point_2(tmp_point->first),tmp_point->second));
          else
            *res++=vec.at(i);
        }

        return res;
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

	std::vector<Object> vec;
	
	CK().intersect_2_object()(c1.arc(),c2.arc(),std::back_inserter(vec)); 

 	for(unsigned i=0; i<vec.size() ; i++)
	{
	  const Rcirc_arc_2 *tmp_arc;

	  if ( (tmp_arc=object_cast<Rcirc_arc_2>(&vec.at(i)) )!=NULL )
	    *res++ = make_object( Circular_arc_2(*tmp_arc) );
	  else
	    {
              const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;
              tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(&vec.at(i));
              CGAL_assertion(tmp_point!=NULL);
              *res++ = make_object( std::make_pair(Circular_arc_point_2(tmp_point->first),tmp_point->second));   
	    }

        }

	return res;
	

      }

     template < class OutputIterator >
    OutputIterator
    operator()(const Line_arc_2 & c1, const Line_arc_2 & c2, 
	       OutputIterator res)
      { 
         Bbox_2 bb1=c1.bbox(),bb2=c2.bbox();

        if(!do_overlap(bb1,bb2 ))
         return res;

	std::vector<Object> vec;
	
	CK().intersect_2_object()(c1.arc(),c2.arc(),std::back_inserter(vec)); 

 	for(unsigned i=0; i<vec.size() ; i++)
	{
	  const Rline_arc_2 *tmp_arc;

	  if ( (tmp_arc=object_cast<Rline_arc_2>(&vec.at(i)) )!=NULL )
	    *res++ = make_object( Line_arc_2(*tmp_arc) );
	  else
            {
              const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;
              tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(&vec.at(i));
              CGAL_assertion(tmp_point!=NULL);
	      *res++ = make_object( std::make_pair(Circular_arc_point_2(tmp_point->first),tmp_point->second));
            }
	}
	return res;	
      }

     template < class OutputIterator >
    OutputIterator
    operator()(const Circular_arc_2 & c1, const Line_arc_2 & c2, 
	       OutputIterator res)
      { 
         Bbox_2 bb1=c1.bbox(),bb2=c2.bbox();

        if(!do_overlap(bb1,bb2 ))
         return res;

	std::vector<Object> vec;
	
	CK().intersect_2_object()(c1.arc(),c2.arc(),std::back_inserter(vec)); 

 	for(unsigned i=0; i<vec.size() ; i++)
	{
	  const Rcirc_arc_2 *tmp_arc;
	  const Rline_arc_2 *tmp_line;

	  if ( (tmp_arc=object_cast<Rcirc_arc_2>(&vec.at(i)) )!=NULL )  // Can this happen?
	    *res++ = make_object( Circular_arc_2(*tmp_arc) );
	  else if ( (tmp_line=object_cast<Rline_arc_2>(&vec.at(i)) )!=NULL ) //Can this happen?
	    *res++ = make_object( Line_arc_2(*tmp_line) );
	  else
	    {
              const std::pair<Rcirc_arc_point_2, unsigned> *tmp_point;
              tmp_point=object_cast<std::pair<Rcirc_arc_point_2, unsigned> >(&vec.at(i));
              CGAL_assertion(tmp_point!=NULL);
              *res++ = make_object( std::make_pair(Circular_arc_point_2(tmp_point->first),tmp_point->second));
            }

	}

	return res;	

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


  template < class BK >
  class Split_2 : public BK::Circular_kernel::Split_2
  {

    typedef typename BK::Circular_kernel       CK;
    typedef typename CK::Circular_arc_2        Rcirc_arc_2;
    typedef typename CK::Line_arc_2            Rline_arc_2;
    typedef typename BK::Circular_arc_2        Circular_arc_2;
    typedef typename BK::Line_arc_2            Line_arc_2;
    typedef typename BK::Circular_arc_point_2  Circular_arc_point_2;

  public:
    typedef typename CK::Split_2::result_type result_type; 
    using CK::Split_2::operator();

    result_type
    operator()(const Circular_arc_2 &A, 
	       const Circular_arc_point_2 &p,
	       Circular_arc_2 &ha1, Circular_arc_2 &ha2) const
    {  
      Rcirc_arc_2 ca1 , ca2;

      CK().split_2_object()(A.arc(), p.point(), ca1, ca2);

      ha1=Circular_arc_2(ca1); 
      ha2=Circular_arc_2(ca2);
    }
    
    result_type
    operator()(const Line_arc_2 &A, 
	       const Circular_arc_point_2 &p,
	       Line_arc_2 &ha1, Line_arc_2 &ha2) const
    {  
      Rline_arc_2 ca1 , ca2;

      CK().split_2_object()(A.arc(), p.point(), ca1, ca2);

      ha1=Line_arc_2(ca1); 
      ha2=Line_arc_2(ca2);
    }
    
  };

template <class BK>
class Construct_bbox_2 : public BK::Circular_kernel::Construct_bbox_2
  {
	  typedef typename BK::Circular_kernel           CK;
    typedef typename BK::Circular_arc_2            Circular_arc_2;
    typedef typename BK::Circular_arc_point_2      Circular_arc_point_2;
    typedef typename BK::Line_arc_2                Line_arc_2;

  public:

	  typedef typename CK::Construct_bbox_2::result_type result_type; 
    using CK::Construct_bbox_2::operator();

    result_type operator() (const Circular_arc_point_2 & a) const
    {
      return a.rep().bbox();
    }

    result_type operator() (const Circular_arc_2 & a) const
    {
      return a.rep().bbox();
    }

    result_type operator() (const Line_arc_2 & a) const
    {
      return a.rep().bbox();
    }

  };


  
} //Bbox_functors

CGAL_END_NAMESPACE  

#endif // CGAL_BBOX_FILTERED_PREDICATES_H  
