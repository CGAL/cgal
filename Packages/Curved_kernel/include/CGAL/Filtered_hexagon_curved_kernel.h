#ifndef CGAL_FILTERED_HEXAGON_CURVED_KERNEL_H
#define CGAL_FILTERED_HEXAGON_CURVED_KERNEL_H


#include <CGAL/Filtered_hexagon_curved_kernel/Circular_arc_with_hexagon_2.h>
#include <CGAL/Filtered_hexagon_curved_kernel/hexagon_filtered_predicates.h>


CGAL_BEGIN_NAMESPACE


template <class CK>
 class Filtered_hexagon_curved_kernel {

  public:

    typedef Filtered_hexagon_curved_kernel<CK>       Self;
    typedef Circular_arc_with_hexagon_2<CK>          Circular_arc_2;
    typedef CK                                       Curved_kernel;
//    typedef typename CK::Linear_kernel               Linear_kernel;
//    typedef typename CK::Algebraic_kernel            Algebraic_kernel;
    typedef typename CK::RT                          RT;
    typedef typename CK::FT                          FT;
    typedef typename CK::Root_of_2                   Root_of_2;
//    typedef typename CK::Polynomial_for_circles_2_2  Polynomial_for_circles_2_2;
//    typedef typename CK::Polynomial_1_2              Polynomial_1_2;
    typedef typename CK::Line_2                      Line_2;
    typedef typename CK::Circle_2                    Circle_2;
    typedef typename CK::Conic_2                     Conic_2;
    typedef typename CK::Point_2                     Point_2;
    typedef typename CK::Circular_arc_2              Rcirc_arc_2;
    typedef typename CK::Circular_arc_point_2        Circular_arc_point_2;
//    typedef typename CK::Line_arc_2                  Line_arc_2;
    typedef typename CK::Construct_circle_2          Construct_circle_2;
    typedef typename CK::Get_equation                Get_equation;

  

    typedef typename CK::Compare_x_2                 Compare_x_2;
    typedef typename CK::Compare_y_2		     Compare_y_2;
    typedef typename CK::Compare_xy_2		     Compare_xy_2;
    typedef CGALi::Construct_Circular_min_vertex_2<Self>   Construct_Circular_min_vertex_2;
    typedef CGALi::Construct_Circular_max_vertex_2<Self>   Construct_Circular_max_vertex_2;
    typedef CGALi::Compare_y_at_x_2<Self>	     Compare_y_at_x_2;
    typedef CGALi::Compare_y_to_right_2<Self>	     Compare_y_to_right_2;
    typedef CGALi::Do_overlap_2<Self>		     Do_overlap_2;
    typedef CGALi::Equal_2<Self>		     Equal_2;
    typedef CGALi::In_range_2<Self>		     In_range_2;
    typedef CGALi::Make_x_monotone_2<Self>	     Make_x_monotone_2;
    typedef CGALi::Intersect_2<Self>                 Intersect_2;
    typedef CGALi::Split_2<Self>		     Split_2;
    typedef CGALi::Is_vertical_2<Self>               Is_vertical_2;






	Get_equation
	get_equation_object() const
	{ return CK().get_equation_object(); }

	Construct_circle_2
	construct_circle_2_object() const
	{ return CK().construct_circle_2_object(); }

	Compare_x_2
	compare_x_2_object() const
	{ return CK().compare_x_2_object(); }

	Compare_y_2
	compare_y_2_object() const
  	{ return CK().compare_y_2_object(); }

	Compare_xy_2
  	compare_xy_2_object() const
    	{ return CK().compare_xy_2_object(); }

	Construct_Circular_min_vertex_2
	construct_circular_min_vertex_2_object() const
  	{ return Construct_Circular_min_vertex_2(); }

	Construct_Circular_max_vertex_2
	construct_circular_max_vertex_2_object() const
  	{ return Construct_Circular_max_vertex_2(); }

  	Compare_y_at_x_2
  	compare_y_at_x_2_object() const 
  	{ return Compare_y_at_x_2(); }

  	Compare_y_to_right_2
  	compare_y_to_right_2_object() const
  	{ return Compare_y_to_right_2(); }

  	Do_overlap_2
  	do_overlap_2_object() const
  	{ return Do_overlap_2(); }

  	Equal_2
  	equal_2_object() const
  	{ return Equal_2(); }

  	In_range_2
  	in_range_2_object() const
  	{ return In_range_2(); }

  	Make_x_monotone_2
  	make_x_monotone_2_object() const
  	{ return Make_x_monotone_2(); }

  	Intersect_2
  	intersect_2_object() const
    	{ return Intersect_2(); }


  	Split_2
  	split_2_object() const
  	{ return Split_2(); }

	Is_vertical_2
	  is_vertical_2_object() const
	{ return Is_vertical_2(); }
};

CGAL_END_NAMESPACE

#endif // CGAL_FILTERED_HEXAGON_CURVED_KERNEL_H
