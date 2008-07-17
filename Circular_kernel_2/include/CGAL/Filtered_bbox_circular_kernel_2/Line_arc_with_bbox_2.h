// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Monique Teillaud, Sylvain Pion, Constantinos Tsirogiannis

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_LINE_ARC_BBOX_HEXAGON_2_H
#define CGAL_LINE_ARC_BBOX_HEXAGON_2_H

#include <vector>
#include <CGAL/Bbox_2.h>
#include <CGAL/Filtered_bbox_circular_kernel_2/Circular_arc_endpoint_with_bbox_2.h>

CGAL_BEGIN_NAMESPACE

template < class BK >
class Line_arc_with_bbox_2 {

    typedef typename BK::Circular_kernel                         CK;
    typedef typename CK::FT                                    FT;
    typedef typename CK::RT                                    RT;
    typedef typename CK::Point_2                               Point_2;
    typedef typename CK::Line_2                                Line_2;
    typedef typename CK::Segment_2                             Segment_2;
    typedef typename CK::Circle_2                              Circle_2;
    typedef typename BK::Circular_arc_point_2                  Circular_arc_point_2;
    typedef Line_arc_with_bbox_2<BK>                           Self;
    typedef typename CK::Line_arc_2                            Line_arc_2;
    typedef typename CK::Root_of_2                             Root_of_2;
    typedef CK R;


public:

                ///////////Construction/////////////

    		Line_arc_with_bbox_2(){}

    		Line_arc_with_bbox_2(const Line_2 &support, 
                       		        const Circle_2 &l1, const bool b_l1,
                       		        const Circle_2 &l2, const bool b_l2)
    		: P_arc(support,l1,b_l1,l2,b_l2), bb(NULL)
    		{}

    
    		Line_arc_with_bbox_2(const Line_2 &support, 
		       		        const Line_2 &l1,
		       		        const Line_2 &l2)
    		: P_arc(support,l1,l2), bb(NULL)
    		{}

		Line_arc_with_bbox_2(const Line_2 &support,
                 		        const Circular_arc_point_2 &begin,
                 		        const Circular_arc_point_2 &end)
    		: P_arc(support, begin.point(), end.point()) , bb(NULL)
		{}


    		Line_arc_with_bbox_2(const Segment_2 &s)
    		: P_arc(s) , bb(NULL)
    		{}

    
    		Line_arc_with_bbox_2(const Point_2 &p1,
                 		     const Point_2 &p2)
    		: P_arc(p1,p2) , bb(NULL)
    		{}

  
		Line_arc_with_bbox_2(const Line_arc_2 &a)
    		: P_arc(a) , bb(NULL)
		{}

    // This avoids Memory Leaks, but may decrease the performance
    // probably not the best solution
    Line_arc_with_bbox_2(const Line_arc_with_bbox_2 &c) : P_arc(c.P_arc), bb(NULL) { }
	  ~Line_arc_with_bbox_2() { if(bb) delete bb; }


		//////////Predicates//////////

		bool is_vertical() const
		{ return P_arc.is_vertical();}

		//////////Accessors///////////

		const Line_arc_2& arc () const
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
		
		const Line_2 & supporting_line() const
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
		
		
	private:

		Line_arc_2 P_arc;
		mutable Bbox_2 *bb;


}; // Line_arc_with_hexagon_2



  template < typename CK >
  std::ostream &
  operator<<(std::ostream & os, const Line_arc_with_bbox_2<CK> &a)
  {
    // The output format is :
    // Supporting line
    // Source endpoint
    // Target endpoint
    return os << a.arc() << " ";
  }

  template < typename CK >
  std::istream &
  operator>>(std::istream & is, Line_arc_with_bbox_2<CK> &a)
  {
    typename CK::Line_2 s;
    typename CK::Circular_arc_point_2 p1;
    typename CK::Circular_arc_point_2 p2;
    is >> s >> p1 >> p2 ;
    if (is)
      a = Line_arc_with_bbox_2<CK>(s, p1, p2);
    return is;
  }

CGAL_END_NAMESPACE

#endif // CGAL_LINE_ARC_WITH_BBOX_2_H
