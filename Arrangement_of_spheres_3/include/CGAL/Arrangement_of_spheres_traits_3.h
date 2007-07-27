#ifndef ARRANGEMENT_OF_SPHERES_TRAITS_3_H
#define ARRANGEMENT_OF_SPHERES_TRAITS_3_H

#include <CGAL/Arrangement_of_spheres_3_basic.h>


#include <CGAL/Root_of_traits.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_3_table.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_line_intersection.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Arrangement_of_spheres_3/Event_point_3.h>
/*#include <CGAL/Arrangement_of_spheres_traits_3/Filtered_sphere_line_intersection.h>
#include <CGAL/Arrangement_of_spheres_traits_3/Sphere_key.h>
#include <CGAL/Arrangement_of_spheres_traits_3/template_types.h>*/

/*#include <CGAL/Arrangement_of_spheres_traits_3/Function_kernel.h>
#include <CGAL/Kinetic/IO/Qt_widget_2.h>
#include <CGAL/Kinetic/Two_list_pointer_event_queue.h>
#include <CGAL/Kinetic/Default_simulator.h>*/

/*#include <CGAL/Arrangement_of_spheres_3/coordinates.h>
  #include <CGAL/Tools/Coordinate_index.h>*/

CGAL_BEGIN_NAMESPACE
#ifdef CGAL_AOS3_USE_TEMPLATES
template <class GT>
#endif
struct Arrangement_of_spheres_traits_3 {
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

  typedef int Event_key;

  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Sphere_line_intersection<This> Sphere_point_3;
  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Event_point_3<This> Event_point_3;
  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Sphere_3_table CGAL_AOS3_TARG Table;
  typedef CGAL_AOS3_TYPENAME Table::Key Sphere_3_key;

  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Coordinate_index Coordinate_index;


  Geom_traits geom_traits_object() const {
    return table_->geom_traits_object();
  }
  
  
  
  
  template <class It> 
  Arrangement_of_spheres_traits_3(It bs, It es): table_(new Table(bs, es)){
    //di_= table_->geometric_traits_object().intersect_3_object();
    std::cout << *table_ << std::endl;
  }
  
  CGAL_GET(Bbox_3, bbox_3, return table_->bbox_3());
  CGAL_SET(Sphere_3, temp_sphere, table_->set_temp_sphere(k));
  CGAL_GET(FT, max_coordinate, return table_->inf());
  CGAL_CONST_ITERATOR(Sphere_3, sphere_3, CGAL_AOS3_TYPENAME Table::Sphere_3_const_iterator,
		      return table_->sphere_3s_begin(),
		      return table_->sphere_3s_end());

  //CGAL_CONST_FIND(Sphere_3, return table_->find(k));
  Sphere_3 sphere_3(Sphere_3_key k) const {
    return table_->operator[](k);
  }
  //CGAL_INSERTNK(Sphere_3, return table_->insert(s));
  CGAL_SIZE(sphere_3s, return table_->number_of_sphere_3s());
 
  /*
    Helpers -----------------------------------------------------------------
  */
  //Plane_3 rule_plane(Key a, Coordinate_index C) const;
  

  /* 
     quadratic constructions ------------------------------------------------
  */

  typedef  std::pair<Event_point_3, Event_point_3> Event_pair;
  
  Event_pair sphere_events(Sphere_3_key s) const;
  /*
  Event_pair intersection_2_events(Key a, Key b) const;

  Event_pair intersection_3_events(Key a, Key b, Key c) const;

  Event_pair sphere_intersect_extremum_events(Key a,  Coordinate_index C,
					      Key b) const;

  Event_pair sphere_intersect_rule_events(Key a, Key r, Coordinate_index C) const;

  Event_pair circle_cross_rule_events(Key a, Key b,
				      Key rs, Coordinate_index C) const;

  Event_pair sphere_intersect_rule_rule_events(Key s, Key rx, Key ry) const;

  // not really used
  Quadratic_NT intersection_c(Key s, Line_3 l, Coordinate_index C) const;*/


  /* 
     predictes--------------------------------------------------------------
  */

  /*bool intersects(Key a, Key b) const;

  bool intersects(Key a, Key b, Key c) const;

  bool sphere_intersects_rule(Key sphere, Key rule_sphere, 
			      Coordinate_index C) const;

  bool sphere_intersects_rule(Key sphere, const Sphere_point_3& sp, 
			      Coordinate_index C) const;
  */
  bool sphere_intersects_sweep(Sphere_3_key sphere, const Sphere_point_3& ep) const;
  /*

  CGAL::Comparison_result compare_equipower_point_to_rule(Key sphere0,
							  Key sphere1,
							  const Sphere_point_3& sp,
							  Coordinate_index C) const;

  CGAL::Sign sign_of_separating_plane_normal_c(Key sphere0, Key sphere1,
					       Coordinate_index C) const ;


  CGAL::Sign sign_of_equipower_plane_normal_c(Key sphere0, Key sphere1, 
					      Coordinate_index C) const;


  CGAL::Oriented_side oriented_side_of_equipower_plane(Key sphere_0, Key sphere_1,
						       const Sphere_point_3 &s) const;
  */
  CGAL::Oriented_side oriented_side_of_separating_plane(Sphere_3_key sphere_0, Sphere_3_key sphere_1, 
							const Sphere_point_3& sp) const ;
  
  CGAL::Comparison_result compare_sphere_centers_c(Sphere_3_key a, Sphere_3_key b, Coordinate_index C) const;

  CGAL::Comparison_result compare_sphere_center_c(Sphere_3_key a,
						  FT d,
						  Coordinate_index C) const;

  CGAL::Comparison_result compare_sphere_center_c(Sphere_3_key a,
						  const Sphere_point_3& d,
						  Coordinate_index C) const;
  // return true if the interval defined by the sphere in C contains d.C
  bool is_over_circle_c(Sphere_3_key c, const Sphere_point_3& d,
			      Coordinate_index C) const;
  /*

  CGAL::Comparison_result compare_sphere_sphere_at_sweep(const Sphere_point_3 &t,
							 Sphere_3_key sphere0,
							 Sphere_3_key Sphere1,
							 const Sphere_point_3 &pt,
							 Coordinate_index C) const;

  // name sucks
  // Find the line which has as coordinate C pt[C] and gets it other coord
  // from planex. See if it is inside the sphere at t
  CGAL::Bounded_side bounded_side_of_sphere_projected( const Sphere_point_3 &t,
						       Sphere_3_key sphere,
						       Sphere_3_key planex,
						       const Sphere_point_3 &pt,
						       Coordinate_index C) const;
  */
  CGAL::Bounded_side bounded_side_of_sphere(Sphere_3_key sphere,
					    const Sphere_point_3 &z) const;
  CGAL::Bounded_side bounded_side_of_sphere(Sphere_3_key sphere,
					    Sphere_3_key x,
					    Sphere_3_key y,
					    const Sphere_point_3 &z) const;
  
  CGAL::Comparison_result compare_depths(const Sphere_point_3 &a, 
					 const Sphere_point_3 &b) const;
  
  CGAL::Oriented_side oriented_side(const Plane_3 &p,
				    const Sphere_point_3 &pt) const;


  // Sweep types ------------------------------------------------------

  /*typedef Arrangement_of_spheres_3_internal::Function_kernel<Event_point_3> Function_kernel;

  typedef CGAL::Kinetic::Default_simulator<Function_kernel,
					   CGAL::Kinetic::Two_list_pointer_event_queue<Function_kernel, true> > Simulator;
					   typedef CGAL::Kinetic::Qt_widget_2<Simulator> Qt_gui;*/
  
private:

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
