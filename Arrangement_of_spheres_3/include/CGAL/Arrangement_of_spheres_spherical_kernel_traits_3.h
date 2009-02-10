#ifndef ARRANGEMENT_OF_SPHERES_SK_TRAITS_3_H
#define ARRANGEMENT_OF_SPHERES_SK_TRAITS_3_H

#include <CGAL/Arrangement_of_spheres_3_basic.h>


#include <CGAL/Arrangement_of_spheres_3/Event_point_3.h>
#include <CGAL/Arrangement_of_spheres_3/Function_kernel.h>
#include <CGAL/Arrangement_of_spheres_3/Rule_direction.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_3_table.h>
#include <CGAL/Arrangement_of_spheres_3/Center_point_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Kinetic/Default_simulator.h>
#include <CGAL/Kinetic/Two_list_pointer_event_queue.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Tools/Label.h>


#define CGAL_AOS3_SK_CANONICAL_PREDICATE_LINE(l, sphere, t, equal_value, inside_value, below_value, above_value, miss_value) \
  std::cout << "Line is " << l << std::endl;				\
  Event_pair ep;							\
  if (!event_pair(sphere, l)) {						\
    std::cout << "Miss " << std::endl;					\
    miss_value;								\
  } else {								\
    Coordinate_index sc=CGAL_AOS3_INTERNAL_NS::sweep_coordinate();	\
    Comparison_result c0= compare_c(t, sp0.first, sc);			\
    Comparison_result c1= compare_c(t, sp1.first, sc);			\
    if (c0 == CGAL::EQUAL || c1 == CGAL::EQUAL) {			\
      std::cout << "Equal" << std::endl;				\
      equal_value;							\
    } else if (c0 != c1) {						\
      std::cout << "Inside" << std::endl;				\
      inside_value;							\
    } else if (c0 == CGAL::SMALLER) {					\
      std::cout << "smaller" << std::endl;				\
      below_value;							\
    } else {								\
      std::cout << "larger" << std::endl;				\
      above_value;							\
    }									\
  }									\
  		
							
#define CGAL_AOS3_SK_CANONICAL_PREDICATE(plane0, plane1, sphere, t, equal_value, inside_value, below_value, above_value, miss_value) \
  Plane_3 p0= plane0;							\
  std::cout << "Plane0 is " << p0 << std::endl;				\
  Plane_3 p1= plane1;							\
  std::cout << "Plane1 is " << p1 << std::endl;				\
  std::vector<Object> out;						\
  i3_(p0,p1, std::back_inserter(out));					\
  CGAL_assertion(out.size()==1);					\
  Line_3 l;								\
  if (assign(l,out.front()) ){						\
    CGAL_AOS3_CANONICAL_PREDICATE_LINE(l, sphere, t, equal_value,	\
				       inside_value, below_value,	\
				       above_value, miss_value);	\
  } else {								\
    std::cout << "Planes miss" << std::endl;				\
    miss_value;								\
  }		

CGAL_BEGIN_NAMESPACE
#ifdef CGAL_AOS3_USE_TEMPLATES
template <class ST>
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
  typedef Arrangement_of_spheres_traits_3<ST> This;
  typedef ST Spherical_traits;
  typedef typename ST::Linear_kernel Linear_traits;
  typedef This Traits_t;
#else
  typedef Arrangement_of_spheres_traits_3 This;
  typedef CGAL_AOS3_INTERNAL_NS::Arrangement_of_spheres_3_geom_traits Linear_traits;
  typedef Exact_spherical_kernel_3<Linear_traits, Algebraic_kernel_for_spheres_2_3<Linear_traits::RT> >  Geom_traits;
#endif
 
  typedef CGAL_AOS3_TYPENAME Spherical_traits::FT FT;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Root_of_2 Root_of_2;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Sphere_3 Sphere_3;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Point_3 Point_3;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Plane_3 Plane_3;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Vector_3 Vector_3;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Segment_3 Segment_3;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Line_3 Line_3;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Line_2 Line_2;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Point_2 Point_2;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Vector_2 Vector_2;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Circle_2 Circle_2;
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Segment_2 Segment_2;
  typedef CGAL_AOS3_INTERNAL_NS::Rule_direction Rule_direction;

  
  typedef CGAL_AOS3_TYPENAME Spherical_traits::Circular_arc_point_3 Sphere_point_3;
  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Event_point_3<This> Event_point_3;
  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Sphere_3_table<This> Table;
  typedef CGAL_AOS3_TYPENAME Table::Key Sphere_3_key;
  typedef CGAL_AOS3_INTERNAL_NS::Function_kernel<This> Function_kernel;
  typedef CGAL::Kinetic::Default_simulator<Function_kernel> Simulator;
  typedef CGAL_AOS3_TYPENAME Simulator::Event_key Event_key;
  

  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Coordinate_index Coordinate_index;

  typedef CGAL_AOS3_INTERNAL_NS::Center_point_3<Sphere_3_key, Event_point_3> Center_point_3;

  
  
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


  Geom_traits geom_traits_object() const {
    return table_->geom_traits_object();
  }
  
  
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

  Plane_3 debug_rule_plane(Sphere_3_key a, Coordinate_index C) const{
    return rule_plane(a,C);
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
  template <class T>
  CGAL::Comparison_result compare_to_circle_circle_c(const T &pt,
						     Sphere_3_key sphere0,
						     Sphere_3_key sphere1,
						     Coordinate_index C) const {
    Plane_3 eqp=table_->equipower_plane(sphere0, sphere1);
    CGAL::Comparison_result cr= compare_to_equipower_line_c(pt, sphere0, sphere1, C);
    //Oriented_side ossp= center_oriented_side_of_separating_plane_c(k, sphere0, sphere1, C);
    Plane_3 sp= table_->separating_plane(sphere0, sphere1);
    CGAL::Comparison_result epqc= compare(sp.orthogonal_vector()[C.index()],0);
    if (cr != epqc ) {
      return cr;
    }
    
    CGAL_AOS3_CANONICAL_PREDICATE(eqp,
				  const_c_plane(pt, C),
				  table_->sphere(sphere0),
				  pt,
				  return CGAL::EQUAL,
				  return Comparison_result(-cr),
				  return cr,
				  return cr,
				  return cr);
  }
  
  //
  template <class T>
  CGAL::Comparison_result compare_to_circle_rule_c(const T &pt,
						   Sphere_3_key sphere,
						   Sphere_3_key rule,
						   Coordinate_index rule_coordinate,
						   bool smaller,
						   Coordinate_index C) const {
    if (C== rule_coordinate) return compare_to_rule_c(pt, rule, rule_coordinate);
    
    Comparison_result crrs= smaller? SMALLER: LARGER; //compare_sphere_centers_c(rule, sphere, C);
    Comparison_result crcs= compare_to_rule_c(pt, sphere, C);
    std::cout << "CPCR " << crrs << " " << crcs << std::endl;
    CGAL_assertion(crrs != CGAL::EQUAL);
    if (crrs != crcs) return crcs;
    
    //Comparison_result cr= compare_point_to_rule_c(pt, rule, C);
    
    //Comparison_result cr= compare_point_to_rule_c(pt, sphere, C);
    CGAL_AOS3_SK_CANONICAL_PREDICATE(const_c_plane(pt, C),
				  rule_plane(rule, rule_coordinate),
				  table_->sphere(sphere),
				  pt,
				  return CGAL::EQUAL,
				  return Comparison_result(-crrs),
				  return crrs,
				  return crrs,
				  return crrs);
  }



  //
 

  template <class T>
  CGAL::Comparison_result compare_to_circle_extremum_c(const T&pt,
						       Sphere_3_key sphere,
						       Rule_direction d,
						       Coordinate_index C) const {
    if (C== d.constant_coordinate()) return compare_to_rule_c(pt, sphere, C);
    
    Comparison_result crrs= d.is_positive()? CGAL::LARGER : CGAL::SMALLER;
    Comparison_result crcs= compare_to_rule_c(pt, sphere, C);
    //CGAL_assertion(crrs != CGAL::EQUAL);
    if (crrs != crcs) return crcs;
    
    //Comparison_result cr= compare_point_to_rule_c(pt, rule, C);
    
    Bounded_side bs= bounded_side_of_sphere_c(pt, sphere, C);
    std::cout << "For " << pt << " and " << sphere << " with RD " << d 
	      << " and direction " << C << " got bs of " << bs << std::endl;
    if (bs== ON_BOUNDARY) return EQUAL;
    if (bs== ON_BOUNDED_SIDE && d.is_positive()
	|| bs == ON_UNBOUNDED_SIDE && d.is_negative()) return SMALLER;
    else return LARGER;
    //Comparison_result cr= compare_point_to_rule_c(pt, sphere, C);
    
  }
  
  
  //
  CGAL::Comparison_result compare_to_rule_c(const Sphere_point_3 &pt,
					    Sphere_3_key rule,
					    Coordinate_index rule_coordinate) const;
  
  CGAL::Comparison_result compare_to_rule_c(const Center_point_3 &pt,
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
  template <class T>
  CGAL::Comparison_result compare_to_equipower_line_c(const T& sp,
						      Sphere_3_key a,
						      Sphere_3_key b,
						      Coordinate_index C) const {
    // NOTE redo this with computing it directly
    CGAL_assertion(C!=CGAL_AOS3_INTERNAL_NS::sweep_coordinate());
    Line_3 l= table_->equipower_line(a,b);
    const Line_3 plt= line_through(sp);
    Line_2 l2(Point_2(l.point()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()],
		      l.point()[C.index()]),
	      Vector_2(l.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()],
		       l.to_vector()[C.index()]));
    Line_2 plt2(Point_2(plt.point()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()],
			plt.point()[C.index()]),
		Vector_2(plt.to_vector()[CGAL_AOS3_INTERNAL_NS::sweep_coordinate().index()],
			 plt.to_vector()[C.index()]));
    //std::cout << "l2 is " << l2 << std::endl;
    //std::cout << "pl2 is " << plt2 << std::endl;
    Object o= geom_traits_object().intersect_2_object()(l2, plt2);
    Point_2 p;
    FT cmp=0;
    Comparison_result c0= compare(plt2.y_at_x(cmp), l2.y_at_x(cmp));
    if (c0 == EQUAL) {
      cmp=1;
      c0=  compare(plt2.y_at_x(cmp), l2.y_at_x(cmp));
      if (c0 == EQUAL) {
	return EQUAL;
      }
    }
    if (assign(p,o)) {
      
      //std::cout << "Sign changes at " << p.x() << std::endl;
      if (compare_c(sp, p.x(),CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) == EQUAL) return EQUAL;
      if (p.x()>0 && compare_c(sp, p.x(), CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) == LARGER
	  || p.x()<0 && compare_c(sp, p.x(), CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) == SMALLER ){
	return Comparison_result(-c0);
      } else return c0;
    } else {
      return c0;
    }
  }
 


  //
  CGAL::Comparison_result compare_c(const Sphere_point_3& a,
				    const Sphere_point_3& b,
				    Coordinate_index C) const {
    return a.compare(b, C);
  }

  CGAL::Comparison_result compare_c(const Center_point_3& a,
				    const Center_point_3& b,
				    Coordinate_index C) const {
    if (C != CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) {
      return compare_sphere_centers_c(a.key(),b.key(),C);
    } else {
      a.coord().compare(b.coord(), C);
    }
  }

  CGAL::Comparison_result compare_c(const Sphere_point_3& a,
				    const Center_point_3& b,
				    Coordinate_index C) const {
    if (C != CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) {
      return compare_to_rule_c(a,b.key(),C);
    } else {
      return a.compare(b.coord(), C);
    }
  }

  CGAL::Comparison_result compare_c(const Center_point_3& a,
				    const Sphere_point_3& b,
				    Coordinate_index C) const {
    if (C != CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) {
      return Comparison_result(-compare_to_rule_c(b,a.key(),C));
    } else {
      return a.coord().compare(b, C);
    }
  }
  
  CGAL::Comparison_result compare_c(const Sphere_point_3& a,
				    const FT& b,
				    Coordinate_index C) const {
    return a.compare(b, C);
  }

  CGAL::Comparison_result compare_c(const Center_point_3& a,
				    const FT& b,
				    Coordinate_index C) const {
    if (C != CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) {
      return CGAL::compare(table_->center(a.key())[C.index()], b);
    } else {
      return a.coord().compare(b, C);
    }
  }
  

  //
  template <class T>
  CGAL::Bounded_side bounded_side_of_sphere(const T &pt,
					    Sphere_3_key s) const {
    
    // This is a bit silly but makes it uniform
    CGAL_AOS3_SK_CANONICAL_PREDICATE(const_c_plane(pt, CGAL_AOS3_INTERNAL_NS::plane_coordinate(0)),
				  const_c_plane(pt, CGAL_AOS3_INTERNAL_NS::plane_coordinate(1)),
				  table_->sphere(s),
				  pt,
				  CGAL_assertion(l== line_through(pt) 
						 || l.opposite() ==line_through(pt));
				  return ON_BOUNDARY,
				  CGAL_assertion(l== line_through(pt) 
						 || l.opposite() ==line_through(pt)); 
				  return ON_BOUNDED_SIDE,
				  CGAL_assertion(l== line_through(pt) 
						 || l.opposite() ==line_through(pt)); 
				  return ON_UNBOUNDED_SIDE,
				  CGAL_assertion(l== line_through(pt) 
						 || l.opposite() ==line_through(pt)); 
				  return ON_UNBOUNDED_SIDE,
				  CGAL_assertion(l== line_through(pt) 
						 || l.opposite() ==line_through(pt)); 
				  return ON_UNBOUNDED_SIDE);
  }
  
  /*CGAL::Bounded_side extremum_bounded_side_of_sphere(const Event_point_3 &t,
    Sphere_3_key pt,
    Rule_direction rd,
    Sphere_3_key sphere) const;*/
  
  //
  CGAL::Bounded_side rules_bounded_side_of_sphere(const Sphere_point_3 &t,
						  Sphere_3_key x,
						  Sphere_3_key y,
						  Sphere_3_key sphere) const;
  
  
  
  template <class T>
  CGAL::Bounded_side bounded_side_of_sphere_c(const T &pt,
					      Sphere_3_key s,
					      Coordinate_index C) const {
    if (C== CGAL_AOS3_INTERNAL_NS::sweep_coordinate()) {
      // also silly but uniform
      
      CGAL_AOS3_SK_CANONICAL_PREDICATE(rule_plane(s, CGAL_AOS3_INTERNAL_NS::plane_coordinate(0)),
				    rule_plane(s, CGAL_AOS3_INTERNAL_NS::plane_coordinate(1)),
				    table_->sphere(s),
				    pt,
				    return ON_BOUNDARY,
				    return ON_BOUNDED_SIDE,
				    return ON_UNBOUNDED_SIDE,
				    return ON_UNBOUNDED_SIDE,
				    return ON_UNBOUNDED_SIDE);
      /*Event_pair ep= sphere_events(s);
	if (pt.compare(ep.first, C)==EQUAL || pt.compare(ep.second, C)==EQUAL) return CGAL::ON_BOUNDARY;
	else if (pt.compare(ep.first, C) == LARGER && pt.compare(ep.second,C)==SMALLER) return CGAL::ON_BOUNDED_SIDE;
	else return CGAL::ON_UNBOUNDED_SIDE;*/
    } else {
      CGAL_AOS3_SK_CANONICAL_PREDICATE(const_c_plane(pt, C),
				    rule_plane(s, CGAL_AOS3_INTERNAL_NS::other_plane_coordinate(C)),
				    table_->sphere(s),
				    pt,
				    return ON_BOUNDARY,
				    return ON_BOUNDED_SIDE,
				    return ON_UNBOUNDED_SIDE,
				    return ON_UNBOUNDED_SIDE,
				    return ON_UNBOUNDED_SIDE);
    }
  }
  
 
  
  CGAL::Comparison_result compare_sphere_centers_c(Sphere_3_key a, Sphere_3_key b, Coordinate_index C) const;

  
  template <class T>
  CGAL::Oriented_side oriented_side_of_separating_plane(const T& sp,
							Sphere_3_key a, 
							Sphere_3_key b) const {
    Plane_3 p= table_->separating_plane(a,b);
    return oriented_side(p, sp);
  }
  
  
 
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
  template <class T>
  CGAL::Oriented_side oriented_side(const Plane_3 &p,
				    const T &s) const {
    CGAL_assertion(0);
  }


  Event_pair event_pair(Sphere_3_key k, const Plane_3 &a, const Plane_3 &b) const;

  Event_pair event_pair(Sphere_3_key k, const Line_3 &l) const ;

  // Sweep types ------------------------------------------------------

  /*typedef Arrangement_of_spheres_3_internal::Function_kernel<Event_point_3> Function_kernel;

    typedef CGAL::Kinetic::Default_simulator<Function_kernel,
    CGAL::Kinetic::Two_list_pointer_event_queue<Function_kernel, true> > Simulator;
    typedef CGAL::Kinetic::Qt_widget_2<Simulator> Qt_gui;*/
  
private:
  


  Event_pair sphere_intersect_rule_rule_events(Sphere_3_key s, Sphere_3_key rx, Sphere_3_key ry) const;

  Event_pair circle_cross_rule_event_internal(Sphere_3_key a, Sphere_3_key b,
					      Sphere_3_key rs, Coordinate_index C) const;

  // plane with C constant
  Plane_3 rule_plane(Sphere_3_key a, Coordinate_index C) const;

  Plane_3 const_c_plane(const Center_point_3 &pt, Coordinate_index C) const {
    return rule_plane(pt.key(), C);
  }


  Plane_3 const_c_plane(const Sphere_point_3 &pt, Coordinate_index C) const;
 
  


  //  HDS hds_;
  // ick, this is to handle location of points which are not already there
  // -1, -2 are bl, tr inf corners
  CGAL_AOS3_TYPENAME Table::Handle table_;
  CGAL_AOS3_TYPENAME Spherical_kernel::Intersect_3 i3_;
  /*mutable Table::Handle table_;
    Geom_traits::Intersect_3 di_;*/
};
CGAL_END_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Arrangement_of_spheres_spherical_kernel_traits_3_impl.h>
#endif
#endif
