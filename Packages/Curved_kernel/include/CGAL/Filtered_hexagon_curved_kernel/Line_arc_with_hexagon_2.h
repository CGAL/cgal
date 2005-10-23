// Author : Constantinos Tsirogiannis

#ifndef CGAL_LINE_ARC_WITH_HEXAGON_2_H
#define CGAL_LINE_ARC_WITH_HEXAGON_2_H

#include <vector>
#include <fstream> //vgale
#include <CGAL/Timer.h> //vgale
#include <iterator>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Curved_kernel/Debug_id.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/NT_extensions_Root_of/CGAL_Interval_nt.h>
#include <CGAL/Filtered_hexagon_curved_kernel/hexagon_primitives.h>

CGAL_BEGIN_NAMESPACE

template < class CK >
class Line_arc_with_hexagon_2 : public CGALi::Debug_id<> {

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
    typedef CGAL::Curved_kernel< FK,CGAL::Algebraic_kernel_2_2<FK::RT> >   CK2;
    typedef CGAL::Curved_kernel_converter<CK,CK2>                          Conv;


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

                typename Qualified_result_of<typename R::Construct_Circular_min_vertex_2,Line_arc_2>::type
		left() const
			{ return P_arc.left();}
                
                typename Qualified_result_of<typename R::Construct_Circular_max_vertex_2,Line_arc_2>::type
                right() const
			{ return P_arc.right();}

		/* const Circular_arc_point_2 & source() const */
/* 			{ return P_arc.source();} */

/* 		const Circular_arc_point_2 & target() const */
/* 			{ return P_arc.target();} */

typename Qualified_result_of<typename R::Construct_Circular_source_vertex_2,Line_arc_2>::type
                source() const
                        {
			  return typename R::Construct_Circular_source_vertex_2()(this->arc());
			}
	      
    typename Qualified_result_of<typename R::Construct_Circular_source_vertex_2,Line_arc_2>::type
                target() const
                        {
			  return typename R::Construct_Circular_target_vertex_2()(this->arc());
			}
		
		const Line_2 & supporting_line() const
			{ return P_arc.supporting_line();}

		Bbox_2 bbox()
			{ return P_arc.bbox();}
			
			
		///Specific check used for hexagon construction///
		
		bool has_no_hexagons() const
		{ return (hxgn==NULL);}
		
		
		///Hexagon construction///
		
		void construct_hexagons() const
		{

                 typedef typename boost::mpl::if_<boost::is_same<typename CK::Definition_tag, typename CK::Curved_tag>, \
                                                  Hexagon_construction_with_interval_2<CK,Hexagon>, \
                                                  Hexagon_construction_on_lazy_kernel_2<CK,Hexagon> >::type Construct;


		  hxgn = new Hexagon(Construct()(P_arc));
		}


	private:

		Line_arc_2 P_arc;
		mutable Hexagon *hxgn;


}; // Line_arc_with_hexagon_2



  template < typename CK >
  std::ostream &
  operator<<(std::ostream & os, const Line_arc_with_hexagon_2<CK> &a)
  {
    // The output format is :
    // Supporting line
    // Source endpoint
    // Target endpoint
    return os << a.arc() << " ";
  }

  template < typename CK >
  std::istream &
  operator>>(std::istream & is, Line_arc_with_hexagon_2<CK> &a)
  {
    typename CK::Line_2 s;
    typename CK::Circular_arc_point_2 p1;
    typename CK::Circular_arc_point_2 p2;
    is >> s >> p1 >> p2 ;
    if (is)
      a = Line_arc_with_hexagon_2<CK>(s, p1, p2);
    return is;
  }

CGAL_END_NAMESPACE

#endif // CGAL_LINE_ARC_WITH_HEXAGON_2_H
