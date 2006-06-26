#ifndef ARRANGEMENT_OF_SPHERES_TRAITS_3_H
#define ARRANGEMENT_OF_SPHERES_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>

#include <CGAL/Root_of_traits.h>
#include <CGAL/Arrangement_of_spheres_3/Filtered_sphere_line_intersection.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_key.h>
#include <CGAL/Tools/Coordinate_index.h>



struct Arrangement_of_spheres_traits_3 {
  typedef Arrangement_of_spheres_traits_3 This;
  typedef CGAL::Cartesian<CGAL::Gmpq> Geometric_traits;

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
  //typedef Filtered_sphere_line_intersection<This, 2> Sphere_point_3;
  typedef Sphere_key Key;
  typedef ::Coordinate_index Coordinate_index;


  Geometric_traits geometric_traits_object() const {
    return t_;
  }

 
  

  template <class It> 
  Arrangement_of_spheres_traits_3(It bs, It es): has_temp_(false){
    spheres_.reserve(std::distance(bs, es)+3);
    initialize_1();
    spheres_.insert(spheres_.end(), bs, es);
    initialize_2();
  }
  
  void set_temp_sphere(const Sphere_3 &s) const {
    has_temp_=true;
    spheres_[0]=s;
  }

  // really just debugging
  Key new_sphere(const Sphere_3 &s) const {
    spheres_.push_back(s);
    return Key(spheres_.size()-4);
  }

  unsigned int number_of_spheres() const {
    return spheres_.size()-3;
  }

  typedef std::vector<Sphere_3>::const_iterator Sphere_iterator;
  Sphere_iterator spheres_begin() const {
    return spheres_.begin()+3;
  }
  Sphere_iterator spheres_end() const {
    return spheres_.end();
  }

  struct Sphere_key_iterator{
    typedef Key value_type;
    typedef const Key &reference_type;
    typedef const Key* pointer_type;
    typedef size_t difference_type;
    Sphere_key_iterator(){}
    Sphere_key_iterator(int i): k_(i){}
    value_type operator*() const {return k_;}
    pointer_type operator->() const {return &k_;}
    Sphere_key_iterator operator++() {
      k_= Key(k_.input_index()+1);
      return *this;
    }
    Sphere_key_iterator operator++(int) {
      Sphere_key_iterator ret=*this;
      operator++();
      return ret;
    }
    bool operator==(const Sphere_key_iterator &o) const {
      return k_==o.k_;
    }
    bool operator!=(const Sphere_key_iterator &o) const {
      return k_!=o.k_;
    }
    Key k_;
  };

  Sphere_key_iterator sphere_keys_begin() const {
    return Sphere_key_iterator(0);
  }
  Sphere_key_iterator sphere_keys_end() const {
    return Sphere_key_iterator(spheres_.size());
  }


  /* 
     linear constructions----------------------------------------------------
  */

  FT inf() const {
    return center_c(Key(Key::TR), Coordinate_index(1));
  }

  // the point described by the vertex (a,b) should be on the positive side
  Plane_3 separating_plane(Key a, Key b) const ;

  // point from the second to the first
  Plane_3 equipower_plane(Key a, Key b) const ;
  
  // point from the second to the first
  Point_3 equipower_point(Key a, Key b) const ;

  FT center_c(Key ind, Coordinate_index C) const;

  Point_3 center(Key ind) const;

  Sphere_3 sphere(Key ind) const;

  

  /* 
     quadratic constructions ------------------------------------------------
  */

  typedef  std::pair<Sphere_point_3, Sphere_point_3> Event_pair;

  Event_pair sphere_events(Key s) const;

  Event_pair intersection_events(Key a, Key b) const;

  Event_pair intersection_events(Key a, Key b, Key c) const;

  // not really used
  Quadratic_NT intersection_c(Key s, Line_3 l, Coordinate_index C) const;


  /* 
     predictes--------------------------------------------------------------
  */

  bool intersects(Key a, Key b) const;

  bool intersects(Key a, Key b, Key c) const;

  bool sphere_intersects_rule(Key sphere, Key rule_sphere, 
			      Coordinate_index C) const;

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
  void initialize_1();
  void initialize_2();
  FT disc(Key i) const;
 

  //  HDS hds_;
  // ick, this is to handle location of points which are not already there
  // -1, -2 are bl, tr inf corners
  mutable std::vector<Sphere_3> spheres_;
  mutable bool has_temp_;
  Geometric_traits t_;
  Geometric_traits::Intersect_3 di_;
};

#endif
