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
// Author(s)     : Monique Teillaud, Sylvain Pion, Julien Hazebrouck, Pedro Machado

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_LINE_ARC_2_H
#define CGAL_LINE_ARC_2_H

CGAL_BEGIN_NAMESPACE

template <class CircularKernel> 
class Line_arc_2 
  : public CircularKernel::Kernel_base::Line_arc_2
{
  typedef typename CircularKernel::FT                        FT;
  typedef typename CircularKernel::RT                        RT;
  //typedef typename CircularKernel::Linear_kernel::Point_2    Point_2;
  typedef typename CircularKernel::Point_2                   Point_2;
  typedef typename CircularKernel::Line_2                    Line_2;
  typedef typename CircularKernel::Circle_2                  Circle_2;
  typedef typename CircularKernel::Circular_arc_point_2   Circular_arc_point_2;
  typedef typename CircularKernel::Segment_2                 Segment_2;

  typedef typename CircularKernel::Kernel_base::Line_arc_2 RLine_arc_2;
public:
  typedef  RLine_arc_2 Rep;
  typedef  CircularKernel   R; 

 const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

   Line_arc_2()
     : RLine_arc_2(typename R::Construct_line_arc_2()())
   {}

   // Not Documented
   Line_arc_2(const Line_2 &support,
	      const Circle_2 &c1,const bool b1,
	      const Circle_2 &c2,const bool b2)
     : RLine_arc_2(typename R::Construct_line_arc_2()(support, c1, b1, c2, b2))
   {}

   // Not Documented
   Line_arc_2(const Line_2 &support,
	       const Line_2 &l1,
	       const Line_2 &l2)
     : RLine_arc_2(typename R::Construct_line_arc_2()(support, l1, l2))
   {}

   Line_arc_2(const Line_2 &support,
	       const Circular_arc_point_2 &p1,
	       const Circular_arc_point_2 &p2)
     : RLine_arc_2(typename R::Construct_line_arc_2()(support, p1, p2))
   {}

   Line_arc_2(const Segment_2 &s)
     : RLine_arc_2(typename R::Construct_line_arc_2()(s))
   {}
   
    Line_arc_2(const Point_2 &p1,
	       const Point_2 &p2)
      : RLine_arc_2(typename R::Construct_line_arc_2()(p1, p2))
   {}

    Line_arc_2(const RLine_arc_2 &a )
     : RLine_arc_2(a)
   {}

  typename Qualified_result_of
  <typename R::Construct_circular_source_vertex_2,Line_arc_2>::type
    //const Circular_arc_point_2 & 
    source() const
  {
        return typename R::Construct_circular_source_vertex_2()(*this);
  }

  typename Qualified_result_of
  <typename R::Construct_circular_target_vertex_2,Line_arc_2>::type
  //const Circular_arc_point_2 & 
    target() const
  {
        return typename R::Construct_circular_target_vertex_2()(*this);
  }

  typename Qualified_result_of
  <typename R::Construct_circular_min_vertex_2,Line_arc_2>::type
  //const Circular_arc_point_2 & left() const
  left() const
  {
        return typename R::Construct_circular_min_vertex_2()(*this);
  }

  typename Qualified_result_of
  <typename R::Construct_circular_max_vertex_2,Line_arc_2>::type
  //const Circular_arc_point_2 & right() const
  right() const
  {
        return typename R::Construct_circular_max_vertex_2()(*this);
  }

  typename Qualified_result_of
  <typename R::Construct_line_2,Line_arc_2>::type
  //const Line_2 & 
    supporting_line() const
  {
        return typename R::Construct_line_2()(*this);
  }
  

  bool is_vertical() const
  {
      return typename R::Is_vertical_2()(*this);
  }
    
  Bbox_2  bbox() const
  {
        return typename R::Construct_bbox_2()(*this);
  }

 };

template < typename CircularKernel >
inline
bool
operator==(const Line_arc_2<CircularKernel> &p,
           const Line_arc_2<CircularKernel> &q)
{
  return CircularKernel().equal_2_object()(p, q);
}

template < typename CircularKernel >
inline
bool
operator!=(const Line_arc_2<CircularKernel> &p,
           const Line_arc_2<CircularKernel> &q)
{
  return ! (p == q);
}


 template < typename CK >
    std::ostream &
    operator<<(std::ostream & os, const Line_arc_2<CK> &a)
    {
      
      return os << a.supporting_line() << " "
		<< a.source() << " "
		<< a.target() << " ";
    }

  template < typename CK >
  std::istream &
  operator>>(std::istream & is, Line_arc_2<CK> &a)
  {
    typename CK::Line_2 l;
    typename CK::Circular_arc_point_2 p1;
    typename CK::Circular_arc_point_2 p2;
    is >> l >> p1 >> p2 ;
    if (is)
      a = Line_arc_2<CK>(l, p1, p2);
    return is;
  }

template < class CK >
struct Filtered_bbox_circular_kernel_2;

template < typename CK >
class Line_arc_2 < Filtered_bbox_circular_kernel_2 < CK > > {

	  typedef Filtered_bbox_circular_kernel_2 < CK >         BK;
    typedef Line_arc_2< BK >                               Self;
    typedef typename CK::FT                                FT;
    typedef typename CK::RT                                RT;
    typedef typename CK::Point_2                           Point_2;
    typedef typename CK::Line_2                            Line_2;
    typedef typename CK::Segment_2                         Segment_2;
    typedef typename CK::Circle_2                          Circle_2;
    typedef typename BK::Circular_arc_point_2            Circular_arc_point_2;
    typedef typename CK::Line_arc_2                        Rline_arc_2;
    typedef typename CK::Root_of_2                         Root_of_2;
    typedef CK R;


public:

                ///////////Construction/////////////

    		Line_arc_2(){}

    		Line_arc_2(const Line_2 &support, 
                       		        const Circle_2 &l1, const bool b_l1,
                       		        const Circle_2 &l2, const bool b_l2)
    		: P_arc(support,l1,b_l1,l2,b_l2), bb(NULL)
    		{}

    
    		Line_arc_2(const Line_2 &support, 
		       		        const Line_2 &l1,
		       		        const Line_2 &l2)
    		: P_arc(support,l1,l2), bb(NULL)
    		{}

		Line_arc_2(const Line_2 &support,
                 		        const Circular_arc_point_2 &begin,
                 		        const Circular_arc_point_2 &end)
    		: P_arc(support, begin.point(), end.point()) , bb(NULL)
		{}


    		Line_arc_2(const Segment_2 &s)
    		: P_arc(s) , bb(NULL)
    		{}

    
    		Line_arc_2(const Point_2 &p1,
                 		     const Point_2 &p2)
    		: P_arc(p1,p2) , bb(NULL)
    		{}

  
		Line_arc_2(const Rline_arc_2 &a)
    		: P_arc(a) , bb(NULL)
		{}

	  Line_arc_2(const Line_arc_2 &c) : P_arc(c.P_arc)
	  {
		  if(c.bb) bb = new Bbox_2(*(c.bb));
			else bb = NULL;	
		}
		
	  ~Line_arc_2() { if(bb) delete bb; }


		//////////Predicates//////////

		bool is_vertical() const
		{ return P_arc.is_vertical();}

		//////////Accessors///////////

		const Rline_arc_2& arc () const
			{ return P_arc ;}
  
		///Interface of the inner arc/// 

                typename Qualified_result_of<typename BK::Construct_circular_min_vertex_2,Self>::type
		left() const
	                {return typename BK::Construct_circular_min_vertex_2()(*this);}
                
                typename Qualified_result_of<typename BK::Construct_circular_max_vertex_2,Self>::type
                right() const
	                {return typename BK::Construct_circular_max_vertex_2()(*this);}

                typename Qualified_result_of<typename BK::Construct_circular_source_vertex_2,Self>::type
                source() const
                        {return typename BK::Construct_circular_source_vertex_2()(*this);}
	      
                typename Qualified_result_of<typename BK::Construct_circular_target_vertex_2,Self>::type
                target() const
                        {return typename BK::Construct_circular_target_vertex_2()(*this);}
		
		const Line_2& supporting_line() const
			{ return P_arc.supporting_line();}

                Bbox_2 bbox() const
                        {
                          if(bb==NULL)
                            bb=new Bbox_2(P_arc.bbox());

                          return *bb;
                        }

			
		///Specific check used for bbox construction///
		
		bool has_no_bbox() const
		{ return (bb==NULL);}
		
	bool equal_ref(const Line_arc_2 &c) const
  {
    return CGAL::identical(P_arc, c.P_arc);      
  }
		
	private:

		Rline_arc_2 P_arc;
		mutable Bbox_2 *bb;


};

CGAL_END_NAMESPACE

#endif // CGAL_LINE_ARC_2_H

   
