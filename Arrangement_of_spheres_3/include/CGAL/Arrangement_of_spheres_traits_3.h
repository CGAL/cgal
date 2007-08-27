#ifndef ARRANGEMENT_OF_SPHERES_TRAITS_3_H
#define ARRANGEMENT_OF_SPHERES_TRAITS_3_H

#include <CGAL/Arrangement_of_spheres_3_basic.h>


#include <CGAL/Root_of_traits.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_3_table.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_line_intersection.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Arrangement_of_spheres_3/Event_point_3.h>
#include <CGAL/Tools/Label.h>
/*#include <CGAL/Arrangement_of_spheres_traits_3/Filtered_sphere_line_intersection.h>
#include <CGAL/Arrangement_of_spheres_traits_3/Sphere_key.h>
#include <CGAL/Arrangement_of_spheres_traits_3/template_types.h>*/

#include <CGAL/Kinetic/Two_list_pointer_event_queue.h>
#include <CGAL/Kinetic/Default_simulator.h>
#include <CGAL/Arrangement_of_spheres_3/Function_kernel.h>

/*#include <CGAL/Arrangement_of_spheres_3/coordinates.h>
  #include <CGAL/Kinetic/IO/Qt_widget_2.h>
  #include <CGAL/Tools/Coordinate_index.h>*/

CGAL_BEGIN_NAMESPACE
#ifdef CGAL_AOS3_USE_TEMPLATES
template <class GT>
#endif
struct Arrangement_of_spheres_traits_3 {
  struct Degeneracy_exception {
    
    template <class T>
    void new_face(T){}

    template <class T>
    void new_edge(T){}

    template <class T>
    void new_vertex(T){}
  };

#ifdef CGAL_AOS3_USE_TEMPLATES
  typedef Arrangement_of_spheres_traits_3<GT> This;
  typedef GT Geom_traits;
  typedef This Traits_t;
#else
  typedef Arrangement_of_spheres_traits_3 This;
  typedef CGAL_AOS3_INTERNAL_NS::Arrangement_of_spheres_3_geom_traits GT;
  typedef GT Geom_traits;
#endif
 
  typedef CGAL_AOS3_TYPENAME Geom_traits::FT FT;
  typedef CGAL_AOS3_TYPENAME CGAL::Root_of_traits<FT>::RootOf_2 Quadratic_NT;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Sphere_3 Sphere_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Point_3 Point_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Plane_3 Plane_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Vector_3 Vector_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Segment_3 Segment_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Line_3 Line_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Line_2 Line_2;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Point_2 Point_2;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Vector_2 Vector_2;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Circle_2 Circle_2;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Segment_2 Segment_2;

  
  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Sphere_line_intersection<This> Sphere_point_3;
  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Event_point_3<This> Event_point_3;
  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Sphere_3_table<This> Table;
  typedef CGAL_AOS3_TYPENAME Table::Key Sphere_3_key;
  typedef CGAL_AOS3_INTERNAL_NS::Function_kernel<This> Function_kernel;
  typedef CGAL::Kinetic::Default_simulator<Function_kernel> Simulator;
  typedef CGAL_AOS3_TYPENAME Simulator::Event_key Event_key;
  

  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Coordinate_index Coordinate_index;


  Geom_traits geom_traits_object() const {
    return table_->geom_traits_object();
  }
  
  
  
  
  template <class It> 
  Arrangement_of_spheres_traits_3(It bs, It es): table_(new Table(bs, es)){
    //di_= table_->geometric_traits_object().intersect_3_object();
    //std::cout << *table_ << std::endl;
  }

  Arrangement_of_spheres_traits_3(){}
  
  CGAL_GET(Bbox_3, bbox_3, return table_->bbox_3());
  CGAL_SET(Sphere_3, temp_sphere, table_->set_temp_sphere(k));
  CGAL_GET(FT, max_coordinate, return table_->inf());
  CGAL_CONST_ITERATOR(Sphere_3, sphere_3, CGAL_AOS3_TYPENAME Table::Sphere_3_const_iterator,
		      return table_->sphere_3s_begin(),
		      return table_->sphere_3s_end());

  CGAL_CONST_ITERATOR(Sphere_3_key, sphere_3_key, CGAL_AOS3_TYPENAME Table::Sphere_key_const_iterator, 
		      return table_->sphere_keys_begin(),
		      return table_->sphere_keys_end());


  //CGAL_CONST_FIND(Sphere_3, return table_->find(k));
  Sphere_3 sphere_3(Sphere_3_key k) const {
    return table_->operator[](k);
  }

  Sphere_3_key new_sphere_3(const Sphere_3 &s);

  //CGAL_INSERTNK(Sphere_3, return table_->insert(s));
  CGAL_SIZE(sphere_3s, return table_->number_of_sphere_3s());

  Point_3 min_corner() const {
    return Point_3(table_->bbox_3().xmin(),
		   table_->bbox_3().ymin(),
		   table_->bbox_3().zmin());
  }
  Point_3 max_corner() const {
    return Point_3(table_->bbox_3().xmax(),
		   table_->bbox_3().ymax(),
		   table_->bbox_3().zmax());
  }
 

  Plane_3 debug_separating_plane(Sphere_3_key a, Sphere_3_key b) const {
    return table_->separating_plane(a,b);
  }

  Plane_3 debug_equipower_plane(Sphere_3_key a, Sphere_3_key b) const {
    return table_->equipower_plane(a,b);
  }

  /*
    Helpers -----------------------------------------------------------------
  */
  

  /* 
     quadratic constructions ------------------------------------------------
  */

  typedef  std::pair<Event_point_3, Event_point_3> Event_pair;
  
  Event_pair sphere_events(Sphere_3_key s) const;
  
  Event_pair intersection_2_events(Sphere_3_key a, Sphere_3_key b) const;
  
  Event_pair intersection_3_events(Sphere_3_key a, Sphere_3_key b, Sphere_3_key c) const;
  
  Event_pair sphere_intersect_extremum_events(Sphere_3_key a,  Coordinate_index C,
					      Sphere_3_key b) const;
  /*
  Event_pair sphere_intersect_rule_events(Sphere_3_key a, Sphere_3_key r, Coordinate_index C) const;
  */

  
  void advance_circle_cross_rule_event(Sphere_3_key a, Sphere_3_key b,
				       Sphere_3_key rs, Coordinate_index C);

  Event_point_3 circle_cross_rule_event(Sphere_3_key a, Sphere_3_key b,
					Sphere_3_key rs, Coordinate_index C) const;

  void advance_sphere_intersect_rule_rule_event(Sphere_3_key s, Sphere_3_key rx, Sphere_3_key ry);

  Event_point_3 sphere_intersect_rule_rule_event(Sphere_3_key s, Sphere_3_key rx, Sphere_3_key ry) const;

 
  // not really used
  /*Quadratic_NT intersection_c(Sphere_3_key s, Line_3 l, Coordinate_index C) const;*/


  /* 
     predictes--------------------------------------------------------------
  */

  bool intersects(Sphere_3_key a, Sphere_3_key b) const;

  bool intersects(Sphere_3_key a, Sphere_3_key b, Sphere_3_key c) const;
  
  /*bool sphere_intersects_rule(Sphere_3_key sphere, Sphere_3_key rule_sphere, 
    Coordinate_index C) const;*/
  
  /*
    returns true if the point is the the slab defined by the top and bottom (in C) 
    of the sphere intersecting the sweep plane
  */
  /*bool 
  compare_point_to_circle_c(const Event_point_3 &point,
			    Sphere_3_key sphere_for_slab,
			    Coordinate_index C) const;

  // warning, they switched
  bool compare_center_to_circle_c(Sphere_3_key sphere_for_point, 
				  const Event_point_3 &t,
				  const Sphere_3_key sphere_for_sphere,
				  Coordinate_index C) const;*/
  

  // 
  CGAL::Comparison_result compare_point_to_circle_circle_c(const Sphere_point_3 &pt,
						     Sphere_3_key sphere0,
						     Sphere_3_key sphere1,
						     Coordinate_index C) const;
  
  //
  CGAL::Comparison_result compare_center_to_circle_circle_c(Sphere_3_key k,
							    const Event_point_3 &t,
							    Sphere_3_key sphere0,
							    Sphere_3_key sphere1,
							    Coordinate_index C) const;
  //
  CGAL::Comparison_result compare_point_to_circle_rule_c(const Sphere_point_3 &pt,
							 Sphere_3_key sphere,
							 Sphere_3_key rule,
							 Coordinate_index rule_coordinate,
							 Coordinate_index C) const;
  //
  CGAL::Comparison_result compare_center_to_circle_rule_c(Sphere_3_key center,
							  const Sphere_point_3 &pt,
							  Sphere_3_key sphere,
							  Sphere_3_key rule,
							  Coordinate_index rule_coordinate,
							  Coordinate_index C) const;
  
  //
  CGAL::Comparison_result compare_point_to_rule_c(const Sphere_point_3 &pt,
						  Sphere_3_key rule,
						  Coordinate_index rule_coordinate) const;
  

 
  //bool sphere_intersects_sweep(Sphere_3_key sphere, const Sphere_point_3& ep) const;
  

  /* CGAL::Comparison_result compare_point_to_equipower_point_c(const Sphere_point_3& sp,
							     Sphere_3_key sphere0,
							     Sphere_3_key sphere1,
							     Coordinate_index C) const;

  CGAL::Comparison_result compare_center_to_equipower_point_c(Sphere_3_key cen,
							      Sphere_3_key sphere0,
							      Sphere_3_key sphere1,
							      Coordinate_index C) const;*/
  //
  CGAL::Comparison_result compare_point_to_equipower_line_c(const Sphere_point_3& sp,
							    Sphere_3_key sphere0,
							    Sphere_3_key sphere1,
							    Coordinate_index C) const;
  //
  CGAL::Comparison_result compare_center_to_equipower_line_c(Sphere_3_key cen,
							     const Event_point_3 &t,
							     Sphere_3_key sphere0,
							     Sphere_3_key sphere1,
							     Coordinate_index C) const;


  //
  CGAL::Comparison_result compare_points_c(const Sphere_point_3& a,
					   const Sphere_point_3& b,
					   Coordinate_index C) const;

  //
  CGAL::Bounded_side point_bounded_side_of_sphere(const Sphere_point_3 &pt,
						  Sphere_3_key sphere) const;
  
  //
  CGAL::Bounded_side center_bounded_side_of_sphere(const Event_point_3 &t,
						   Sphere_3_key pt,
						   Sphere_3_key sphere) const;
  
  //
  CGAL::Bounded_side rules_bounded_side_of_sphere(const Sphere_point_3 &t,
						  Sphere_3_key x,
						  Sphere_3_key y,
						  Sphere_3_key sphere) const;


  
  CGAL::Bounded_side point_bounded_side_of_sphere_c(const Sphere_point_3 &pt,
						    Sphere_3_key sphere,
						    Coordinate_index C) const;
  
  CGAL::Bounded_side center_bounded_side_of_circle_c(Sphere_3_key k,
						    const Event_point_3 &t,
						    Sphere_3_key sphere,
						    Coordinate_index C) const;


 
  
  CGAL::Comparison_result compare_sphere_centers_c(Sphere_3_key a, Sphere_3_key b, Coordinate_index C) const;

  
  CGAL::Oriented_side point_oriented_side_of_separating_plane(const Sphere_point_3& sp,
							      Sphere_3_key sphere_0, 
							      Sphere_3_key sphere_1) const ;
 
  CGAL::Oriented_side center_oriented_side_of_separating_plane(Sphere_3_key center,
							       const Sphere_point_3& t,
							       Sphere_3_key sphere_0, 
							       Sphere_3_key sphere_1) const ;
  

protected:

 
protected:
  /*
  bool sphere_over_rule(Sphere_3_key sphere, const Sphere_point_3& sp, 
			      Coordinate_index C) const;
 
  
 

  CGAL::Comparison_result compare_sphere_center_c(Sphere_3_key a,
						  FT d,
						  Coordinate_index C) const;*/

  /*CGAL::Comparison_result compare_sphere_center_c(Sphere_3_key a,
						  const Sphere_point_3& d,
						  Coordinate_index C) const;*/


  // return true if the interval defined by the sphere in C contains d.C
   
 

  /*CGAL::Comparison_result compare_sphere_sphere_at_sweep(const Sphere_point_3 &t,
							 Sphere_3_key sphere0,
							 Sphere_3_key Sphere1,
							 const Sphere_point_3 &pt,
							 Coordinate_index C) const;*/
  
 

  // name sucks
  // Find the line which has as coordinate C pt[C] and gets it other coord
  // from planex. See if it is inside the sphere at t
  /*CGAL::Bounded_side bounded_side_of_sphere_projected( const Sphere_point_3 &t,
						       Sphere_3_key sphere,
						       Sphere_3_key planex,
						       const Sphere_point_3 &pt,
						       Coordinate_index C) const;
 
  CGAL::Bounded_side bounded_side_of_sphere(Sphere_3_key sphere,
					    const Sphere_point_3 &z) const;
  CGAL::Bounded_side bounded_side_of_sphere(Sphere_3_key sphere,
					    Sphere_3_key x,
					    Sphere_3_key y,
					    const Sphere_point_3 &z) const;
  
  CGAL::Comparison_result compare_depths(const Sphere_point_3 &a, 
					 const Sphere_point_3 &b) const;
  */
  CGAL::Oriented_side oriented_side(const Plane_3 &p,
				    const Sphere_point_3 &pt) const;


  // Sweep types ------------------------------------------------------

  /*typedef Arrangement_of_spheres_3_internal::Function_kernel<Event_point_3> Function_kernel;

  typedef CGAL::Kinetic::Default_simulator<Function_kernel,
					   CGAL::Kinetic::Two_list_pointer_event_queue<Function_kernel, true> > Simulator;
					   typedef CGAL::Kinetic::Qt_widget_2<Simulator> Qt_gui;*/
  
private:
  Event_pair sphere_intersect_rule_rule_events(Sphere_3_key s, Sphere_3_key rx, Sphere_3_key ry) const;

  Event_pair circle_cross_rule_events(Sphere_3_key a, Sphere_3_key b,
				      Sphere_3_key rs, Coordinate_index C) const;

  Plane_3 rule_plane(Sphere_3_key a, Coordinate_index C) const;

  Plane_3 const_c_plane(const Sphere_point_3 &pt, Coordinate_index C) const;
 

  //  HDS hds_;
  // ick, this is to handle location of points which are not already there
  // -1, -2 are bl, tr inf corners
  CGAL_AOS3_TYPENAME Table::Handle table_;

  /*mutable Table::Handle table_;
    Geom_traits::Intersect_3 di_;*/
};
CGAL_END_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Arrangement_of_spheres_traits_3_impl.h>
#endif
#endif
