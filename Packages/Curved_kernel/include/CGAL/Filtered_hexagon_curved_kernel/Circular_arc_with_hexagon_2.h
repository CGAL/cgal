// Author : Constantinos Tsirogiannis

#ifndef CGAL_CIRCULAR_ARC_WITH_HEXAGON_2_H
#define CGAL_CIRCULAR_ARC_WITH_HEXAGON_2_H

#include <vector>
#include <iterator>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Curved_kernel/Debug_id.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Filtered_hexagon_curved_kernel/hexagon_primitives.h>

CGAL_BEGIN_NAMESPACE

template < class CK, 
	   class Container = std::vector<Polygon_2<Simple_cartesian<double> > > >
class Circular_arc_with_hexagon_2 : public CGALi::Debug_id<> {

    typedef typename CK::FT                                    FT;
    typedef typename CK::RT                                    RT;
    typedef typename CK::Linear_kernel::Point_2                Point_2;
    typedef typename CK::Line_2                                Line_2;
    typedef typename CK::Circle_2                              Circle_2;
    typedef typename CK::Circular_arc_point_2               Circular_arc_point_2;
    typedef typename CK::Circular_arc_2                        Circular_arc_2;
    typedef typename CK::Root_of_2                             Root_of_2;

public:

                ///////////TYPE DEFINITIONS/////////////

    typedef typename Container::iterator                       Hexagon_iterator;	
    typedef typename Container::const_iterator                 Hexagon_const_iterator;	
    typedef Polygon_2< Simple_cartesian<double> >              Hexagon;	


                ///////////Construction/////////////

    		Circular_arc_with_hexagon_2(){}

    		Circular_arc_with_hexagon_2(const Circle_2 &c)
    		: P_arc(c)
    		{}

    		Circular_arc_with_hexagon_2(const Circle_2 &support, 
                       		   const Line_2 &l1, const bool b_l1,
                       		   const Line_2 &l2, const bool b_l2)
    		: P_arc(support,l1,b_l1,l2,b_l2)
    		{}

    
    		Circular_arc_with_hexagon_2(const Circle_2 &c, 
		       		   const Circle_2 &c1, const bool b_1,
		       		   const Circle_2 &c2, const bool b_2)
    		: P_arc(c,c1,b_1,c2,b_2)
    		{}

    
    		Circular_arc_with_hexagon_2(const Circular_arc_2 &A, const bool b,
		       		   const Circle_2 &ccut, const bool b_cut)
    		: P_arc(A, b, ccut, b_cut)
    		{}


    		Circular_arc_with_hexagon_2(const Point_2 &start,
                 		   const Point_2 &middle,
                 		   const Point_2 &end)
    		: P_arc(start, middle, end) 
    		{}

  
    		Circular_arc_with_hexagon_2(const Circle_2 &support,
                 		   const Point_2 &begin,
                 		   const Point_2 &end)
    		: P_arc(support, begin, end) 
    		{}


		Circular_arc_with_hexagon_2(const Circle_2 &support,
                 		   const Circular_arc_point_2 &begin,
                 		   const Circular_arc_point_2 &end)
    		: P_arc(support, begin, end) 
		{}

		Circular_arc_with_hexagon_2(const Circular_arc_2 &a)
    		: P_arc(a) 
		{}



		//////////Predicates//////////

		bool is_x_monotone() const
		{ return P_arc.is_x_monotone();}

		bool is_y_monotone() const
		{ return P_arc.is_y_monotone();}

		bool on_upper_part() const
		{ return P_arc.on_upper_part();}
		
		
		//////////Accessors///////////

		const Circular_arc_2& arc () const
			{ return P_arc ;}
  
		Hexagon_iterator hexagons_begin()  
			{ return hexagons.begin();} 

		Hexagon_iterator hexagons_end() 
			{ return hexagons.end();} 
              

		Hexagon_const_iterator hexagons_begin() const 
			{ return hexagons.begin();} 

		Hexagon_const_iterator hexagons_end() const
			{ return hexagons.end();} 

		unsigned hexagon_no() const
			{ return hexagons.size();} 

		
		///Interface of the inner arc/// 

		const Circular_arc_point_2 & left() const
			{ return P_arc.left();}

		const Circular_arc_point_2 & right() const
			{ return P_arc.right();}

		const Circular_arc_point_2 & source() const
			{ return P_arc.source();}

		const Circular_arc_point_2 & target() const
			{ return P_arc.target();}
		
		const Circle_2 & supporting_circle() const
			{ return P_arc.supporting_circle();}

		const Point_2 & center() const
			{ return P_arc.center();}

		const FT & squared_radius() const
			{ return P_arc.squared_radius();}
		
		Bbox_2 bbox()
			{ return P_arc.bbox();}
			
			
		///Specific check used for hexagon construction///
		
		bool has_no_hexagons() const
		{ return hexagons.empty();}
		
		
		///Hexagon construction///
		
		void construct_hexagons() 
		{
		  assert(has_no_hexagons());	
		  CGALi::construct_bounding_hexagons_2<CK>(P_arc,std::back_inserter(hexagons));
		}


	private:

		Circular_arc_2 P_arc;
		Container hexagons;


}; // Circular_arc_with_hexagon_2



  template < typename CK >
  std::ostream &
  operator<<(std::ostream & os, const Circular_arc_with_hexagon_2<CK> &a)
  {
    // The output format is :
    // - supporting circle
    // - circle c1
    // - bool b1
    // - circle c2
    // - bool b2
    return os << a.arc() << " ";
  }

  template < typename CK >
  std::istream &
  operator>>(std::istream & is, Circular_arc_with_hexagon_2<CK> &a)
  {
    typename CK::Circle_2 s;
    typename CK::Circular_arc_point_2 p1;
    typename CK::Circular_arc_point_2 p2;
    is >> s >> p1 >> p2 ;
    if (is)
      a = Circular_arc_with_hexagon_2<CK>(s, p1, p2);
    return is;
  }

  template < typename CK >
  std::ostream &
  print(std::ostream & os, const Circular_arc_with_hexagon_2<CK> &a)
  {
    return os << "Circular_arc_2( " << a.id() << std::endl
              << "left : " << a.arc().left() << " , " << std::endl
              << "right : " << a.arc().right() << " , " << std::endl
	      << "upper part : " << a.arc().on_upper_part() << std::endl
              << "  [[ approximate circle is (x,y,r) : "
              << to_double(a.arc().supporting_circle().center().x()) << "  "
              << to_double(a.arc().supporting_circle().center().y()) << "  "
              << std::sqrt(to_double(a.arc().supporting_circle()
                                            .squared_radius()))
              << " ]]" << std::endl;
  }

CGAL_END_NAMESPACE

#endif // CGAL_CIRCULAR_ARC_WITH_HEXAGON_2_H
