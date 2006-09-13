#ifndef ARRANGEMENT_OF_SPHERES_TRAITS_3_H
#define ARRANGEMENT_OF_SPHERES_TRAITS_3_H

#include <CGAL/basic.h>


#include <CGAL/Root_of_traits.h>
#include <CGAL/Arrangement_of_spheres_3/Filtered_sphere_line_intersection.h>
#include <CGAL/Arrangement_of_spheres_3/Event_point_3.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_key.h>
#include <CGAL/Tools/Coordinate_index.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_3_table.h>

#include <CGAL/Arrangement_of_spheres_3/coordinates.h>


struct Arrangement_of_spheres_traits_3 {
  typedef Arrangement_of_spheres_traits_3 This;
 
  typedef ::Geometric_traits Geometric_traits;
  typedef Geometric_traits::FT FT;
  typedef CGAL::Root_of_traits<FT>::RootOf_2 Quadratic_NT;
  typedef Geometric_traits::Sphere_3 Sphere_3;
  typedef Geometric_traits::Point_3 Point_3;
  typedef Geometric_traits::Plane_3 Plane_3;
  typedef Geometric_traits::Vector_3 Vector_3;
  typedef Geometric_traits::Segment_3 Segment_3;
  typedef Geometric_traits::Line_3 Line_3;
  typedef Geometric_traits::Line_2 Line_2;
  typedef Geometric_traits::Point_2 Point_2;
  typedef Geometric_traits::Vector_2 Vector_2;
  typedef Geometric_traits::Circle_2 Circle_2;
  typedef Geometric_traits::Segment_2 Segment_2;
  typedef Sphere_line_intersection<This> Sphere_point_3;
  typedef Event_point_3<This> Event_point_3;
  //typedef Filtered_sphere_line_intersection<This, 2> Sphere_point_3;
  typedef Sphere_3_table Table;
  typedef Table::Key Key;

  typedef ::Coordinate_index Coordinate_index;


  Geometric_traits geometric_traits_object() const {
    return table_->geometric_traits_object();
  }
  
  
  
  
  template <class It> 
  Arrangement_of_spheres_traits_3(It bs, It es): table_(new Table(bs, es)){
    di_= table_->geometric_traits_object().intersect_3_object();
  }
  
  CGAL::Bbox_3 bbox_3() const {
    return table_->bbox_3();
  }
  
  
  void set_temp_sphere(const Sphere_3 &s) const {
    table_->set_temp_sphere(s);
  }

  FT max_coordinate() const {
    return table_->inf();
  }

  Sphere_3 sphere(Key k) const {
    return table_->sphere(k);
  }

  // really just debugging
  Key new_sphere(const Sphere_3 &s) const {
    return table_->new_sphere(s);
  }

  unsigned int number_of_spheres() const {
    return table_->size();
  }

  typedef Table::Sphere_iterator Sphere_iterator;
  Sphere_iterator spheres_begin() const {
    return table_->spheres_begin();
  }
  Sphere_iterator spheres_end() const {
    return table_->spheres_end();
  }

  typedef Table::Sphere_key_iterator Sphere_key_iterator;

  Sphere_key_iterator sphere_keys_begin() const {
    return table_->sphere_keys_begin();
  }
  Sphere_key_iterator sphere_keys_end() const {
    return table_->sphere_keys_end();
  }


  /*
    Helpers -----------------------------------------------------------------
  */
  Plane_3 rule_plane(Key a, Coordinate_index C) const;
  

  /* 
     quadratic constructions ------------------------------------------------
  */

  typedef  std::pair<Event_point_3, Event_point_3> Event_pair;

  Event_pair sphere_events(Key s) const;

  Event_pair intersection_2_events(Key a, Key b) const;

  Event_pair intersection_3_events(Key a, Key b, Key c) const;

  Event_pair sphere_intersect_extremum_events(Key a,  Coordinate_index C,
					      Key b) const;

  Event_pair sphere_intersect_rule_events(Key a, Key r, Coordinate_index C) const;

  Event_pair circle_cross_rule_events(Key a, Key b,
				      Key rs, Coordinate_index C) const;

  Event_pair sphere_intersect_rule_rule_events(Key s, Key rx, Key ry) const;

  // not really used
  Quadratic_NT intersection_c(Key s, Line_3 l, Coordinate_index C) const;


  /* 
     predictes--------------------------------------------------------------
  */

  bool intersects(Key a, Key b) const;

  bool intersects(Key a, Key b, Key c) const;

  bool sphere_intersects_rule(Key sphere, Key rule_sphere, 
			      Coordinate_index C) const;

  bool sphere_intersects_sweep(Key sphere, Sphere_point_3 ep) const;

  CGAL::Comparison_result compare_equipower_point_to_rule(Key sphere0,
							  Key sphere1,
							  Key rule_sphere,
							  Coordinate_index C) const;

  CGAL::Sign sign_of_separating_plane_normal_c(Key sphere0, Key sphere1,
					       Coordinate_index C) const ;


  CGAL::Sign sign_of_equipower_plane_normal_c(Key sphere0, Key sphere1, 
					      Coordinate_index C) const;


  CGAL::Oriented_side oriented_side_of_equipower_plane(Key sphere_0, Key sphere_1,
						       const Sphere_point_3 &s) const;

  CGAL::Oriented_side oriented_side_of_center_plane(Key sphere_0, Key sphere_1, 
						    Key sphere_center) const ;
 
  CGAL::Comparison_result compare_sphere_centers_c(Key a, Key b, Coordinate_index C) const;

  CGAL::Comparison_result compare_sphere_center_c(Key a,
						  FT d,
						  Coordinate_index C) const;

  CGAL::Comparison_result compare_sphere_center_c(Key a,
						  const Sphere_point_3& d,
						  Coordinate_index C) const;

  CGAL::Bounded_side bounded_side_of_sphere_on_equipower_plane_rule_line(Key sphere0,
									 Key Sphere1,
									 Key rule,
									 Coordinate_index C,
									 const Sphere_point_3 &sp) const;

  CGAL::Bounded_side bounded_side_of_sphere(Key sphere,
					    Key sphere_x,
					    Key sphere_y,
					    const Sphere_point_3 &z) const;

  CGAL::Comparison_result compare_depths(const Sphere_point_3 &a, 
					 const Sphere_point_3 &b) const;
  
private:

  //  HDS hds_;
  // ick, this is to handle location of points which are not already there
  // -1, -2 are bl, tr inf corners

  mutable Table::Handle table_;
  Geometric_traits::Intersect_3 di_;
};

#endif
