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

#ifndef CGAL_LINE_ARC_WITH_HEXAGON_2_H
#define CGAL_LINE_ARC_WITH_HEXAGON_2_H

#include <vector>
#include <iterator>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Filtered_hexagon_circular_kernel_2/hexagon_primitives.h>

namespace CGAL {

template < class HK >
class Line_arc_with_hexagon_2 {

    typedef typename HK::Circular_kernel                         CK;
    typedef typename CK::FT                                    FT;
    typedef typename CK::RT                                    RT;
    typedef typename CK::Point_2                               Point_2;
    typedef typename CK::Line_2                                Line_2;
    typedef typename CK::Segment_2                             Segment_2;
    typedef typename CK::Circle_2                              Circle_2;
    typedef typename CK::Circular_arc_point_2                  Circular_arc_point_2;
    typedef typename CK::Line_arc_2                            Line_arc_2;
    typedef typename CK::Root_of_2                             Root_of_2;
    typedef CK R;
    typedef CGAL::Simple_cartesian<CGAL::Interval_nt<> >                   FK;
    typedef CGAL::Circular_kernel_2< FK,CGAL::Algebraic_kernel_for_circles_2_2<FK::RT> >   CK2;
    typedef CGAL::Circular_kernel_converter<CK,CK2>                          Conv;


public:

                ///////////TYPE DEFINITIONS/////////////

    typedef Line_arc_2                                         Encapsulated_arc_type;
    typedef Polygon_2< Simple_cartesian<double> >              Hexagon;	
    typedef Hexagon*                                           Hexagon_iterator;	
    typedef const Hexagon *                                    Hexagon_const_iterator;	


                ///////////Construction/////////////

    		Line_arc_with_hexagon_2(){}

    		Line_arc_with_hexagon_2(const Line_2 &support, 
                       		        const Circle_2 &l1, const bool b_l1,
                       		        const Circle_2 &l2, const bool b_l2)
    		: P_arc(support,l1,b_l1,l2,b_l2), hxgn(NULL)
    		{}

    
    		Line_arc_with_hexagon_2(const Line_2 &support, 
		       		        const Line_2 &l1,
		       		        const Line_2 &l2)
    		: P_arc(support,l1,l2), hxgn(NULL)
    		{}

		Line_arc_with_hexagon_2(const Line_2 &support,
                 		        const Circular_arc_point_2 &begin,
                 		        const Circular_arc_point_2 &end)
    		: P_arc(support, begin, end) , hxgn(NULL)
		{}


    		Line_arc_with_hexagon_2(const Segment_2 &s)
    		: P_arc(s) , hxgn(NULL)
    		{}

    
    		Line_arc_with_hexagon_2(const Point_2 &p1,
                 		        const Point_2 &p2)
    		: P_arc(p1,p2) , hxgn(NULL)
    		{}

  
		Line_arc_with_hexagon_2(const Line_arc_2 &a)
    		: P_arc(a) , hxgn(NULL)
		{}



		//////////Predicates//////////

		bool is_vertical() const
		{ return P_arc.is_vertical();}

		//////////Accessors///////////

		const Line_arc_2& arc () const
			{ return P_arc ;}
  
		Hexagon_iterator hexagons_begin()  
			{ return hxgn;} 

		Hexagon_iterator hexagons_end() 
			{ return ( (has_no_hexagons())? NULL: hxgn+1);} 

		Hexagon_const_iterator hexagons_begin() const 
			{ return hxgn;} 

		Hexagon_const_iterator hexagons_end() const
			{ return ( (has_no_hexagons())? NULL: hxgn+1);} 

		unsigned hexagon_no() const
			{ return ( (has_no_hexagons())? 0 : 1);} 

		
		///Interface of the inner arc/// 

                typename Qualified_result_of<typename R::Construct_circular_min_vertex_2,Line_arc_2>::type
		left() const
			{ return P_arc.left();}
                
                typename Qualified_result_of<typename R::Construct_circular_max_vertex_2,Line_arc_2>::type
                right() const
			{ return P_arc.right();}

                typename Qualified_result_of<typename R::Construct_circular_source_vertex_2,Line_arc_2>::type
                source() const
                        {
			  return typename R::Construct_circular_source_vertex_2()(this->arc());
			}
	      
                typename Qualified_result_of<typename R::Construct_circular_target_vertex_2,Line_arc_2>::type
                target() const
                        {
			  return typename R::Construct_circular_target_vertex_2()(this->arc());
			}
		
		const Line_2 & supporting_line() const
			{ return P_arc.supporting_line();}

		Bbox_2 bbox() const
			{ return P_arc.bbox();}
			
			
		///Specific check used for hexagon construction///
		
		bool has_no_hexagons() const
		{ return (hxgn==NULL);}
		
		
		///Hexagon construction///
		
		void construct_hexagons() const
		{

                 typedef typename boost::mpl::if_<boost::is_same<typename CK::Definition_tag, typename CK::Circular_tag>, \
                                                  Hexagon_construction_with_interval_2<CK,Hexagon>, \
                                                  Hexagon_construction_on_lazy_kernel_2<CK,Hexagon> >::type Construct;
		  hxgn = new Hexagon(Construct()(P_arc));
		}


	private:

		Line_arc_2 P_arc;
		mutable Hexagon *hxgn;


}; // Line_arc_with_hexagon_2



  template < typename HK >
  std::ostream &
  operator<<(std::ostream & os, const Line_arc_with_hexagon_2<HK> &a)
  {
    // The output format is :
    // Supporting line
    // Source endpoint
    // Target endpoint
    return os << a.arc() << " ";
  }

  template < typename HK >
  std::istream &
  operator>>(std::istream & is, Line_arc_with_hexagon_2<HK> &a)
  {
    typedef typename HK::Circular_kernel                         CK;
    typename CK::Line_2 s;
    typename CK::Circular_arc_point_2 p1;
    typename CK::Circular_arc_point_2 p2;
    is >> s >> p1 >> p2 ;
    if (is)
      a = Line_arc_with_hexagon_2<HK>(s, p1, p2);
    return is;
  }

} //namespace CGAL

#endif // CGAL_LINE_ARC_WITH_HEXAGON_2_H
