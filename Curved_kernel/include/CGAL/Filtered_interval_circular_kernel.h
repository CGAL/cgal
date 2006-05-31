#ifndef FILTERED_INTERVAL_CIRCULAR_KERNEL_H_
#define FILTERED_INTERVAL_CIRCULAR_KERNEL_H_

#include <CGAL/Filtered_interval_circular_kernel/ah_filtered_predicates.h>

CGAL_BEGIN_NAMESPACE

template <class CK>
 class Filtered_interval_circular_kernel : public CK {

  public:

    typedef Filtered_interval_circular_kernel<CK>    Self;
    typedef CK                                       Circular_kernel;
    typedef typename CK::Circular_arc_2              Circular_arc_2;
    typedef typename CK::Line_arc_2                  Line_arc_2;
    typedef typename CK::Circular_arc_point_2        Circular_arc_point_2;
    typedef typename CK::Algebraic_kernel            Algebraic_kernel;
    typedef typename CK::RT                          RT;
    typedef typename CK::FT                          FT;
    typedef typename CK::Root_of_2                   Root_of_2;
    typedef typename CK::Root_for_circles_2_2        Root_for_circles_2_2;
    typedef typename CK::Line_2                      Line_2;
    typedef typename CK::Circle_2                    Circle_2;
    typedef typename CK::Conic_2                     Conic_2;
    typedef typename CK::Point_2                     Point_2;
    typedef typename CK::Circular_arc_2              Rcirc_arc_2;
    typedef typename CK::Line_arc_2                  Rline_arc_2;
    typedef typename CK::Circular_arc_point_2        Rcirc_arc_point_2;
    typedef typename CK::Construct_circle_2          Construct_circle_2;
    typedef typename CK::Get_equation                Get_equation;
    typedef typename CK::Compute_Circular_x_2        Compute_Circular_x_2;
    typedef typename CK::Compute_Circular_y_2        Compute_Circular_y_2;

    typedef AH_functors::Construct_circular_source_vertex_2<Self> 
      Construct_circular_source_vertex_2;
    typedef AH_functors::Construct_circular_target_vertex_2<Self> 
      Construct_circular_target_vertex_2;
    typedef AH_functors::Construct_circular_min_vertex_2<Self> 
      Construct_circular_min_vertex_2;
    typedef AH_functors::Construct_circular_max_vertex_2<Self> 
      Construct_circular_max_vertex_2;
    typedef AH_functors::Compare_x_2<Self> Compare_x_2;
    typedef AH_functors::Compare_y_2<Self> Compare_y_2;
    typedef AH_functors::Compare_xy_2<Self> Compare_xy_2;
    typedef AH_functors::Has_on_2<Self> Has_on_2;
    typedef AH_functors::Compare_y_at_x_2<Self>	Compare_y_at_x_2;
    typedef AH_functors::Compare_y_to_right_2<Self> Compare_y_to_right_2;
    typedef AH_functors::Do_overlap_2<Self> Do_overlap_2;
    typedef AH_functors::Equal_2<Self> Equal_2;
    typedef AH_functors::In_x_range_2<Self> In_x_range_2;
    typedef AH_functors::Make_x_monotone_2<Self> Make_x_monotone_2;
    typedef AH_functors::Intersect_2<Self> Intersect_2;
    typedef AH_functors::Split_2<Self> Split_2;
    typedef AH_functors::Is_vertical_2<Self> Is_vertical_2;

	Get_equation get_equation_object() const { 
	  return CK().get_equation_object(); 
	}

	Construct_circle_2 construct_circle_2_object() const { 
	  return CK().construct_circle_2_object(); 
	}

	Compare_x_2	compare_x_2_object() const { 
	  return Compare_x_2(); 
	}

	Compare_y_2	compare_y_2_object() const { 
	  return Compare_y_2(); 
	}

	Compare_xy_2 compare_xy_2_object() const { 
      return Compare_xy_2(); 
    }

	Construct_circular_source_vertex_2
	construct_circular_source_vertex_2_object() const { 
	  return Construct_circular_source_vertex_2(); 
	}

	Construct_circular_target_vertex_2
	construct_circular_target_vertex_2_object() const { 
	  return Construct_circular_target_vertex_2(); 
	}

	Construct_circular_min_vertex_2
	construct_circular_min_vertex_2_object() const { 
	  return Construct_circular_min_vertex_2(); 
	}

	Construct_circular_max_vertex_2
	construct_circular_max_vertex_2_object() const { 
      return Construct_circular_max_vertex_2(); 
    }

    Has_on_2 has_on_2_object() const { 
      return Has_on_2(); 
    }

  	Compare_y_at_x_2 compare_y_at_x_2_object() const {
  	  return Compare_y_at_x_2(); 
  	}

  	Compare_y_to_right_2
  	compare_y_to_right_2_object() const { 
      return Compare_y_to_right_2(); 
    }

  	Do_overlap_2 do_overlap_2_object() const { 
      return Do_overlap_2(); 
    }

  	Equal_2 equal_2_object() const {
  	  return Equal_2(); 
  	}

  	In_x_range_2 in_x_range_2_object() const { 
  	  return In_x_range_2(); 
  	}

  	Make_x_monotone_2 make_x_monotone_2_object() const { 
  	  return Make_x_monotone_2(); 
  	}

  	Intersect_2 intersect_2_object() const { 
  	  return Intersect_2(); 
  	}

  	Split_2 split_2_object() const { 
      return Split_2(); 
    }

	Is_vertical_2 is_vertical_2_object() const { 
      return Is_vertical_2(); 
    }
};

CGAL_END_NAMESPACE

#endif /*FILTERED_INTERVAL_CIRCULAR_KERNEL_H_*/
