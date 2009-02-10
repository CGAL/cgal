#ifndef CGAL_AOS3_IRRATIONAL_CROSS_SECTION_RULES_H
#define CGAL_AOS3_IRRATIONAL_CROSS_SECTION_RULES_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>
#include <CGAL/Arrangement_of_spheres_3/Rational_cross_section.h>

#define DPRINT(x) 
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Irrational_cross_section_rules {
  CGAL_AOS3_TRAITS;
  typedef Irrational_cross_section_rules CGAL_AOS3_TARG This;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
  //typedef Cross_section_events CGAL_AOS3_TARG CSE;
  //typedef CGAL_AOS3_TYPENAME CS::Halfedge_handle Halfedge_handle;
  //typedef CGAL_AOS3_TYPENAME CS::Face_handle Face_handle;
public:

  typedef CGAL_AOS3_TYPENAME Traits::Sphere_3_key Sphere_3_key;
  typedef CGAL_AOS3_TYPENAME Traits::Event_point_3 Event_point_3;
  typedef CGAL_AOS3_TYPENAME Traits::Sphere_point_3 Sphere_point_3;
  typedef CGAL_AOS3_TYPENAME Traits::Center_point_3 Center_point_3;
  
  struct On_vertex_exception {
    On_vertex_exception(CGAL_AOS3_TYPENAME CS::Vertex_handle h): v_(h){}
    CGAL_AOS3_TYPENAME CS::Vertex_handle vertex_handle() const {
      return v_;
    }
    CGAL_AOS3_TYPENAME CS::Vertex_handle v_;
  };

  Irrational_cross_section_rules(const Traits &tr, CS &cs): tr_(tr),
						      cs_(cs)
						      {}

  


  /*
    If more than one face, make a vertices pass Then make another edges
    pass. For each arc, check if the location of the point differs in
    one direction. If so, compute the intersection circle of the rule
    plane of the point with the sphere. If there is no contact, then the
    edge is OK. If there is contact and the beginning/end of the circle
    surround the insertion point, then we remove the face from
    consideration.
  */



 

 
  CGAL_AOS3_TYPENAME CS::Halfedge_handle 
  find_rule_vertex(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &t, 
		   CGAL_AOS3_TYPENAME CS::Face_handle f,
		   CGAL_AOS3_TYPENAME CS::Curve rule) ;


private:
  
 // return comparison of point on edge of face to the shot rule on the C coordinate
  // i.e. SMALLER if the point on edge is smaller than the rule coordinate
  template <class T>
  CGAL::Comparison_result compare_to_edge_vertex(const T &pt,
						 CGAL_AOS3_TYPENAME CS::Rule_direction rd,
						 CGAL_AOS3_TYPENAME CS::Point p) const {
    bool exact;
    CGAL::Comparison_result debug_answer= debug_compare_to_answer(pt, rd,
								  p, 
								  exact);
    
    Comparison_result cr;
    
    if (compare_to_rational(pt, rd, p, cr)){
    } else if (p.is_sphere_rule()) {
      cr= compare_to_SR(pt, rd, p);
    } else {
      cr= compare_to_SS(pt, rd,p);
    }
    debug_check(debug_answer, cr, exact);
    return cr;
}



  template <class T>
  CGAL::Comparison_result compare_to_SS(const T &t,
					CGAL_AOS3_TYPENAME CS::Rule_direction rd,
					const CGAL_AOS3_TYPENAME CS::Point &point) const  {
    CGAL_AOS3_TYPENAME Traits::Coordinate_index C= rd.constant_coordinate();
    
    CGAL::Comparison_result cr= tr_.compare_to_circle_circle_c(t,
							       point.sphere_key(0),
							       point.sphere_key(1),
							       C);
    return cr;
  }
  


  template <class T>
  CGAL::Comparison_result compare_to_SR(const T &t,
					CGAL_AOS3_TYPENAME CS::Rule_direction rd,
					const CGAL_AOS3_TYPENAME CS::Point &point) const  {
    /*CGAL_LOG(Log::LOTS,"SR " << point << " " << center << " at " << t 
      << ": " << rd << std::endl);*/
    //CGAL_precondition(srule.is_vertical() != orule.is_vertical());
    //CGAL_precondition(arc.is_arc());
    //CGAL_precondition(orule.is_rule());
    //CGAL_precondition(srule.constant_coordinate() != orule.constant_coordinate());
    
    CGAL::Comparison_result cr;
    if (point.is_sphere_extremum()) {
      cr=  tr_.compare_to_circle_extremum_c(t,
					    point.sphere_key(),
					    point.sphere_extremum_index(),
					    rd.constant_coordinate());
    } else {
      cr= tr_.compare_to_circle_rule_c(
				       t,
				       point.sphere_key(),
				       point.rule_key(),
				       point.rule_constant_coordinate(),
				       point.is_smaller(),
				       rd.constant_coordinate());
    }
    
    return cr;
  }




  template <class T>
 bool compare_to_rational(const T & t,
			  CGAL_AOS3_TYPENAME CS::Rule_direction rd,
			  CGAL_AOS3_TYPENAME CS::Point point,
			  CGAL::Comparison_result &ret) const  {
    //NT coord;
    if (point.is_rule_rule()) {
      ret= tr_.compare_to_rule_c(t, point.rule_key(rd.constant_coordinate()),
				 rd.constant_coordinate());
      return true;
    } else if (point.is_sphere_rule() 
	       && point.rule_constant_coordinate() == rd.constant_coordinate()) {
      ret= tr_.compare_to_rule_c(t, point.rule_key(),
				 rd.constant_coordinate());
      return true;
    } else {
      return false;
    }
  }

 
  
  
public:


  template <class T>
  CGAL_AOS3_TYPENAME CS::Halfedge_handle shoot_rule(const T& pt,
						    CGAL_AOS3_TYPENAME CS::Face_handle f,
						    CGAL_AOS3_TYPENAME CS::Rule_direction rd) {
    bool backwards= rd.is_backwards();
    CGAL_LOG(Log::LOTS, "Rule is " << rd << std::endl);
    
    CGAL_AOS3_TYPENAME CS::Halfedge_handle end= f->halfedge();
    
    while (!rd.can_intersect(end->curve())) {
      end=end->next();
    }
    while (rd.can_intersect(end->curve())) {
      end=end->next();
    }
    CGAL_AOS3_TYPENAME CS::Halfedge_handle h=end;
    while (!rd.can_intersect(h->curve())) {
      h=h->next();
    }
    
    do {
      // later skip second one if first is wrong
      CGAL_LOG_WRITE(Log::LOTS, cs_.write(h, LOG_STREAM) 
		     << " is being tested" << std::endl);
      if (h->next() == end) {
	CGAL_LOG_WRITE(Log::LOTS, cs_.write(h, LOG_STREAM) 
		       << " is returned by default" << std::endl);
	return h;
      } else {
	CGAL::Comparison_result cr= compare_to_edge_vertex(pt,
							   rd,
							   h->vertex()->point());
	CGAL_LOG(Log::LOTS, "Result is " << cr << std::endl);
	if (cr == CGAL::EQUAL) {
	  CGAL_LOG(Log::LOTS, "On vertex "
		   << h->vertex()->point() << std::endl);
	  throw On_vertex_exception(h->vertex());
	}
	if (backwards && cr == CGAL::LARGER
	    || !backwards && cr == CGAL::SMALLER) {
	  CGAL_LOG_WRITE(Log::LOTS, cs_.write(h, LOG_STREAM) 
			 << " is returned " << std::endl);
	  return h;
	}
      }
      h= h->next();
    } while (h != end);
    
    CGAL_assertion(0);
    return CGAL_AOS3_TYPENAME CS::Halfedge_handle();
  }

 
 void 
  roll_back_rule(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &t,
		 CGAL_AOS3_TYPENAME CS::Halfedge_handle cur);


protected:

  CGAL::Comparison_result debug_compare_to_answer(const Sphere_point_3&t,
						  CGAL_AOS3_TYPENAME CS::Rule_direction rd,
						  CGAL_AOS3_TYPENAME CS::Point point,
						  bool &exact) const;

  CGAL::Comparison_result debug_compare_to_answer(const Center_point_3&t,
						  CGAL_AOS3_TYPENAME CS::Rule_direction rd,
						  CGAL_AOS3_TYPENAME CS::Point point,
						  bool &exact) const ;


 
  
  void debug_check(CGAL::Comparison_result check, 
		   CGAL::Comparison_result computed,
		   bool exact) const;


  /*bool Slice::rule_shoot_compare_if_rational_arc(const T::Sphere_point_3 & ep,
    Sds::Curve rule,
    Sds::Curve a,
    CGAL::Comparison_result &ret) const {
    CGAL::Comparison_result comp= t_.compare_sphere_centers_c(a.key(), rule.key(),
    rule.constant_coordinate());
    // NOTE could use cached locations
    if (rule.constant_coordinate()==plane_coordinate(1)){
    if (a.is_top() && comp == CGAL::LARGER) {
    ret= CGAL::LARGER;
    return true;
    } else if (!a.is_top() && comp ==CGAL::SMALLER){
    ret= CGAL::SMALLER;
    return true;
    }
    } else {
    if (a.is_right() && comp == CGAL::LARGER) {
    ret= CGAL::LARGER;
    return true;
    } else if (!a.is_right() && comp ==CGAL::SMALLER){
    ret= CGAL::SMALLER;
    return true;
    }
    }
    return false;
    }*/


 




  //! Take the vertex pointed to by h with the face h->face() and put a rule in it.
  CGAL_AOS3_TYPENAME CS::Halfedge_handle
  extend_rule(const CGAL_AOS3_TYPENAME Traits::Event_point_3 &t,
	      CGAL_AOS3_TYPENAME CS::Halfedge_handle h,
	      CGAL_AOS3_TYPENAME CS::Curve rule);
  
  void clean_up_vertex(const CGAL_AOS3_TYPENAME Traits::Event_point_3 &t,
		       CGAL_AOS3_TYPENAME CS::Vertex_handle vh) ;
 

  Traits tr_;
  CS &cs_;
 
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include "Irrational_cross_section_rules_impl.h"
#endif
#endif
