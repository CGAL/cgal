#ifndef CGAL_AOS3_IRRATIONAL_CROSS_SECTION_H
#define CGAL_AOS3_IRRATIONAL_CROSS_SECTION_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>
#include <CGAL/Arrangement_of_spheres_3/Rational_cross_section.h>


#define DPRINT(x) 
CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Irrational_cross_section {
  CGAL_AOS3_TRAITS;
  typedef Irrational_cross_section CGAL_AOS3_TARG This;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
  //typedef CGAL_AOS3_TYPENAME CS::Halfedge_handle Halfedge_handle;
  //typedef CGAL_AOS3_TYPENAME CS::Face_handle Face_handle;
public:

  struct On_edge_exception {
    On_edge_exception(CGAL_AOS3_TYPENAME CS::Halfedge_handle h): h_(h){}
    CGAL_AOS3_TYPENAME CS::Halfedge_handle halfedge_handle() const {
      return h_;
    }
    CGAL_AOS3_TYPENAME CS::Halfedge_handle h_;
  };
  struct On_vertex_exception {
    On_vertex_exception(CGAL_AOS3_TYPENAME CS::Vertex_handle h): v_(h){}
    CGAL_AOS3_TYPENAME CS::Vertex_handle vertex_handle() const {
      return v_;
    }
    CGAL_AOS3_TYPENAME CS::Vertex_handle v_;
  };

  Irrational_cross_section(const Traits &tr, CS &cs, 
			   const CGAL_AOS3_TYPENAME Traits::Event_point_3 &z): tr_(tr),
									       cs_(cs),
									       z_(z){}

  


  /*
    If more than one face, make a vertices pass Then make another edges
    pass. For each arc, check if the location of the point differs in
    one direction. If so, compute the intersection circle of the rule
    plane of the point with the sphere. If there is no contact, then the
    edge is OK. If there is contact and the beginning/end of the circle
    surround the insertion point, then we remove the face from
    consideration.
  */



 

  /*Slice::Face_handle Slice::locate_point(const T::Sphere_point_3 & ep) {
    t_.set_temp_sphere(ep.sphere());
    return locate_point(ep, T::Key::temp_key());
    }*/



  template <class It>
  CGAL_AOS3_TYPENAME CS::Face_handle 
  locate_point(It b, It e,
	       const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & ep) {

    std::vector<int> locations(tr_.number_of_sphere_3s(), 0);
    std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> faces;
    std::vector<CGAL_AOS3_TYPENAME CS::Curve> edges;
    
    {
      std::cout << "Initial faces are " << std::endl;
      for (It c= b; c != e; ++c){
	cs_.write(*c,  std::cout ) << std::endl;
      }
    }

    for (It fit = b; fit != e; ++fit){
      if (!cs_.is_in_slice(*fit)) continue;
      bool ok=locate_point_check_face(ep, *fit, locations/*, edges*/);
      if (ok) faces.push_back(*fit);
    }
    
    {
      std::cout << "After point check " << std::endl;
      for (unsigned int i=0; i< faces.size(); ++i){
	cs_.write(faces[i], std::cout ) << std::endl;
      }
    }

    CGAL_assertion(!faces.empty());
    if (faces.size() > 1) {
      std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> clean_faces;
      for (unsigned int i=0; i< faces.size(); ++i){
	if (locate_point_check_face_arcs(ep, faces[i], locations)) {
	  clean_faces.push_back(faces[i]);
	} 
      }
      std::swap(faces, clean_faces);
    }
    
    {
      std::cout << "After arcs check " << std::endl;
      for (unsigned int i=0; i< faces.size(); ++i){
	cs_.write(faces[i], std::cout ) << std::endl;
      }
    }



    if (faces.size() > 1) {
      std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> clean_faces;
      for (unsigned int i=0; i< faces.size(); ++i){
	if (locate_point_check_face_vertices(ep, faces[i])) {
	  clean_faces.push_back(faces[i]);
	} 
      }
      std::swap(faces, clean_faces);
    }
    
    {
      std::cout << "After vertices check " << std::endl;
      for (unsigned int i=0; i< faces.size(); ++i){
	cs_.write(faces[i], std::cout ) << std::endl;
      }
    }

    if (faces.size() ==1) {
      return faces[0];
    } else if (faces.size() ==2) {
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h= faces[0]->halfedge();
      do {
	if (h->opposite()->face() == faces[1]) throw On_edge_exception(h);
	h= h->next();
      } while (h != faces[0]->halfedge());
      
      cs_.write(faces[0], std::cerr) << std::endl;
      cs_.write(faces[1], std::cerr) << std::endl;
      CGAL_assertion(0);
    } else {
      CGAL_assertion(!faces.empty());
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h= faces[0]->halfedge();
      do {
	bool ok=true;
	for (unsigned int i=1; i< faces.size(); ++i) {
	  if (cs_.has_vertex(faces[i], h->vertex())) {
	    
	  } else { ok=false;}
	}
	if (ok)   throw On_vertex_exception(h->vertex());
	h= h->next();
      } while (h != faces[0]->halfedge());
      
      
      // must handle degeneracy for intersection
      std::cerr << "Remaining faces are " << std::endl;
      for (unsigned int i=0; i< faces.size(); ++i){
	cs_.write(faces[i], std::cerr ) << std::endl;
      }
      CGAL_assertion(0);
      
    }
  
    CGAL_assertion(0);
    return CGAL_AOS3_TYPENAME CS::Face_handle();
  }


  CGAL_AOS3_TYPENAME CS::Face_handle 
  locate_point(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & ep);






  CGAL_AOS3_TYPENAME CS::Halfedge_handle shoot_rule(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t,
						    CGAL_AOS3_TYPENAME CS::Face_handle f,
						    const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & pt,
						    CGAL_AOS3_TYPENAME CS::Rule_direction rd) {
    CGAL_assertion(t.compare(z_, CGAL_AOS3_INTERNAL_NS::Sweep_coordinate::object())== CGAL::EQUAL);
    CGAL_assertion(pt.compare(z_, CGAL_AOS3_INTERNAL_NS::Sweep_coordinate::object())== CGAL::EQUAL);
    //Sds::Curve rule= Sds::Curve::make_rule(ruledir);
    bool backwards= rd.is_backwards();
    DPRINT(std::cout << "Rule is " << rd << std::endl);
  
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
      DPRINT(sds_.write(h, std::cout) << " is being tested" << std::endl);
      if (h->next() == end) {
	DPRINT(sds_.write(h, std::cout) << " is returned by default" << std::endl);
	return h;
      } else {
	CGAL::Comparison_result cr= rule_shoot_edge_vertex(t,
							   pt, rd,
							   h->curve(),
							   h->vertex()->point(),
							   h->next()->curve());
	DPRINT(std::cout << "Result is " << cr << std::endl);
	if (cr == CGAL::EQUAL) {
	  DPRINT(std::cout << "On vertex " << h->vertex()->point() << std::endl);
	  throw On_vertex_exception(h->vertex());
	}
	if (backwards && cr == CGAL::SMALLER
	    || !backwards && cr == CGAL::LARGER) {
	  DPRINT(sds_.write(h, std::cout) << " is returned " << std::endl);
	  return h;
	}
      }
      h= h->next();
    } while (h != end);

    CGAL_assertion(0);
    return CGAL_AOS3_TYPENAME CS::Halfedge_handle();
  }



  // return comparison of point on edge of face to the shot rule on the C coordinate
  // i.e. SMALLER if the point on edge is smaller than the rule coordinate
  CGAL::Comparison_result rule_shoot_edge_vertex(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t,
						 const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & pt, 
						 CGAL_AOS3_TYPENAME CS::Rule_direction rd,
						 CGAL_AOS3_TYPENAME CS::Curve hp,
						 CGAL_AOS3_TYPENAME CS::Point p,
						 CGAL_AOS3_TYPENAME CS::Curve hn) const {
    //CGAL_assertion(h != b);
    //CGAL_assertion(h->curve() != b->curve());
    CGAL::Comparison_result ret;
    // add filters on sphere centers when arcs are involved.
    if (rule_shoot_compare_if_rational(t, pt, rd, p, ret)){
     
      bool exact;
      CGAL::Comparison_result debug_answer= debug_rule_shoot_answer(t, pt, rd,
								    p, 
								    exact);
      debug_rule_shoot_check(debug_answer, ret, exact);
      return ret;
    } else if (p.is_sphere_rule()) {
      if (hp.is_rule() && hn.is_arc()){
	return rule_shoot_compare_SR(t, pt, rd, hn, hp.key(),
				     p,
				     !rd.is_vertical() && 
				     (hn.quadrant() & CS::Curve::T_BIT) 
				     || rd.is_vertical() && 
				     (hn.quadrant() & CS::Curve::R_BIT));
      } else if (hp.is_arc() && hn.is_rule()){
	return rule_shoot_compare_SR(t, pt, rd, hp, hn.key(), 
				     p,
				     !rd.is_vertical() && 
				     (hp.quadrant() & CS::Curve::T_BIT) 
				     || rd.is_vertical() && 
				     (hp.quadrant() & CS::Curve::R_BIT));
      } else if (hp.key() == hn.key()){
	CGAL_assertion(hp.is_arc() && hn.is_arc());
	if (hp == hn) {
	  //CGAL_assertion(hp==hn);
	  bool arc_top;
	  if (rd.constant_coordinate() == plane_coordinate(0)) {
	    arc_top= hp.is_right();
	  } else {
	    arc_top = hp.is_top();
	  }
	  if (!hp.is_inside()) arc_top= !arc_top;
	  
	  return rule_shoot_compare_SR(t, pt, rd, hp, p.rule_key(),p, arc_top);
	} else {
	  // compare rational should have picked this up
	  CGAL_assertion(0);
	}
      }
    } else {
      bool tangent= (hp.key()== hn.key());
      if (tangent) {
	std::pair<CGAL_AOS3_TYPENAME Traits::Event_point_3, 
	  CGAL_AOS3_TYPENAME Traits::Event_point_3> ep3
	  = tr_.intersection_2_events(p.sphere_key(0),
				      p.sphere_key(1));
	std::cout << "Tangent." << std::endl;
	CGAL::Comparison_result cr= ep3.first.compare(t, sweep_coordinate());
	//debug_rule_shoot_check(debug_answer, cr, exact);
	return cr;
      } else {
	return rule_shoot_compare_SS(t, pt, rd,p);
      }
    }

    CGAL_assertion(0);
    return CGAL::SMALLER;
  }


  CGAL::Comparison_result rule_shoot_compare_SS(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t, 
						const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & pt, 
						CGAL_AOS3_TYPENAME CS::Rule_direction rd,
						CGAL_AOS3_TYPENAME CS::Point point) const {
    bool exact;
    CGAL::Comparison_result debug_answer= debug_rule_shoot_answer(t, pt, rd,
								  point,
								  exact);
    //std::cout << "SS " << arc0 << "--" << pt << "--" << arc1 << std::endl;
   
    CGAL_AOS3_TYPENAME Traits::Coordinate_index C= rd.constant_coordinate();

    CGAL::Comparison_result cr= tr_.compare_sphere_sphere_at_sweep(t,
								   point.sphere_key(0),
								   point.sphere_key(1),
								   pt, C);
    debug_rule_shoot_check(debug_answer, cr, exact);
    return cr;
  }



  CGAL::Comparison_result rule_shoot_compare_SR(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t, 
						const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & pt,
						CGAL_AOS3_TYPENAME CS::Rule_direction rd,
						CGAL_AOS3_TYPENAME CS::Curve arc,
						CGAL_AOS3_TYPENAME Traits::Sphere_3_key orule,
						CGAL_AOS3_TYPENAME CS::Point debug_pt,
						bool arc_above) const {
    bool exact;
    CGAL::Comparison_result debug_answer= debug_rule_shoot_answer(t, pt, rd, 
								  debug_pt,
								  exact);
    DPRINT(std::cout << "SR " << arc << "-" << arc_above << "-" << orule 
	   << ": " << rd << std::endl);
    //CGAL_precondition(srule.is_vertical() != orule.is_vertical());
    CGAL_precondition(arc.is_arc());
    //CGAL_precondition(orule.is_rule());
    //CGAL_precondition(srule.constant_coordinate() != orule.constant_coordinate());

    CGAL::Comparison_result c= tr_.compare_sphere_center_c(arc.key(), pt,
							   rd.constant_coordinate());
    if (rd.is_vertical()) {
      if (c== CGAL::LARGER && arc.is_right()) {
	debug_rule_shoot_check(debug_answer, CGAL::LARGER, exact);
	return CGAL::LARGER;
      } else if (c== CGAL::SMALLER && arc.is_left()) {
	debug_rule_shoot_check(debug_answer, CGAL::SMALLER, exact);
	return CGAL::SMALLER;
      }
    } else {
      if (c== CGAL::SMALLER && arc.is_bottom()) {
	debug_rule_shoot_check(debug_answer, CGAL::SMALLER, exact);
	return CGAL::SMALLER;
      } else if (c== CGAL::LARGER && arc.is_top()) {
	debug_rule_shoot_check(debug_answer, CGAL::LARGER, exact);
	return CGAL::LARGER;
      }
    }

      
    CGAL::Bounded_side intersection_side;
    if (rd.constant_coordinate() == plane_coordinate(0)) {
      intersection_side=tr_.bounded_side_of_sphere_projected(t, arc.key(),
							     orule,
							     pt,
							     rd.constant_coordinate());
    } else {
      intersection_side=tr_.bounded_side_of_sphere_projected(t, arc.key(),
							     orule,
							     pt,
							     rd.constant_coordinate());
    }
  
  

    if (intersection_side== CGAL::ON_BOUNDARY){
      debug_rule_shoot_check(debug_answer, CGAL::EQUAL, exact);
      return CGAL::EQUAL;
    } else {
      int cum=0;
      if (arc.is_inside()) ++cum;
      if (intersection_side== CGAL::ON_BOUNDED_SIDE) ++cum;
      if (arc_above) ++cum;
      
      if (cum%2==1) {
	debug_rule_shoot_check(debug_answer, CGAL::LARGER, exact);
	return CGAL::LARGER;
      } else {
	debug_rule_shoot_check(debug_answer, CGAL::SMALLER, exact);
	return CGAL::SMALLER;
      }
    }
    
  }


  CGAL::Comparison_result debug_rule_shoot_answer(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t, 
						  const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & pt,
						  CGAL_AOS3_TYPENAME CS::Rule_direction rd,
						  CGAL_AOS3_TYPENAME CS::Point point,
						  bool &exact) const {
    //CGAL::Comparison_result answer;
    CGAL_AOS3_TYPENAME Traits::FT z;
    if (t.has_simple_coordinate(sweep_coordinate())) {
      z= t.simple_coordinate(sweep_coordinate());
      exact=true;
    } else {
      z= t.approximate_coordinate(sweep_coordinate());
      exact=false;
    }
    /*DPRINT(std::cout << "Comparing " << t_.center(rule.key()) << " on coord " 
      << rule.constant_coordinate() << " which is " 
      << t_.center(rule.key())[rule.constant_coordinate()] << std::endl);*/

    //std::cout << "Z is " << z << " and exact is " << exact << std::endl;

    Rational_cross_section CGAL_AOS3_TARG rcs(cs_, tr_);
    rcs.set_z(z);

    //Sds::Point pt(p,n);
    const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & sp= rcs.sphere_point(point);
    /*if (!sp.is_valid()) {
      const T::Sphere_point_3 & sp2= sphere_point_rz(pt, z);
      }*/
    CGAL_assertion(sp.is_valid());
    
    DPRINT(std::cout << "The point " << pt << " is " 
	   << CGAL::to_double(sp.exact_coordinate(plane_coordinate(0))) << " " 
	   << CGAL::to_double(sp.exact_coordinate(plane_coordinate(1))) << std::endl);
  
    DPRINT(std::cout << "The center is " 
	   << CGAL::to_double(pt.exact_coordinate(plane_coordinate(0))) << " " 
	   << CGAL::to_double(pt.exact_coordinate(plane_coordinate(1))) << std::endl);
  
    CGAL::Comparison_result cr= sp.compare(pt, rd.constant_coordinate());
   
    DPRINT(std::cout << "Got " << cr << std::endl);
    if (cr== CGAL::EQUAL) {
      /*std::cout << "EQUAL for " << p << " " << n << ": " << sp 
	<<  " vs " << ep << " on " 
	<< rule.constant_coordinate() << std::endl;*/
    }
    return cr;
  }

  
  
  void debug_rule_shoot_check(CGAL::Comparison_result check, 
			      CGAL::Comparison_result computed,
			      bool exact) const {
    if (check!= computed) {
      if (exact){
	CGAL_assertion(check==computed);
      } else {
	std::cerr << "Warning, computed " << computed << " expected " 
		  << check << std::endl;
      }
    }
  }


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


  bool rule_shoot_compare_if_rational(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t, 
				      const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & pt,
				      CGAL_AOS3_TYPENAME CS::Rule_direction rd,
				      CGAL_AOS3_TYPENAME CS::Point point,
				      CGAL::Comparison_result &ret) const {
    //NT coord;
    if (point.is_rule_rule()) {
      ret= tr_.compare_sphere_center_c(point.rule_key(rd.constant_coordinate()), pt,
				       rd.constant_coordinate());
      return true;
    } else if (point.is_sphere_rule() 
	       && point.rule_coordinate() == rd.constant_coordinate()) {
      ret= tr_.compare_sphere_center_c(point.rule_key(), pt,
				       rd.constant_coordinate());
      return true;
    } else {
      return false;
    }
  }





private:

  bool locate_point_check_face(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &z,
			       CGAL_AOS3_TYPENAME CS::Face_const_handle it,
			       std::vector<int> &locations) const ;





  // cache checked arcs
  bool locate_point_check_face_arcs(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &ep,
				    CGAL_AOS3_TYPENAME CS::Face_const_handle f,
				    std::vector<int> &locations) const ;



  bool locate_point_check_face_vertices(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &ep,
					CGAL_AOS3_TYPENAME CS::Face_const_handle it) const ;


  int sphere_location(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3& sp,
		      CGAL_AOS3_TYPENAME Traits::Sphere_3_key s) const ;




 




  bool behind_arc(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &ep,
		  CGAL_AOS3_TYPENAME CS::Curve arc,
		  int location) const;


  void point_sphere_orientation(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &time,
				CGAL_AOS3_TYPENAME Traits::Sphere_3_key sphere,
				std::vector<int> &locations
				/*,
				  std::vector<CS::Curve> &edges*/) const ;
  



  Traits tr_;
  CS &cs_;
  CGAL_AOS3_TYPENAME Traits::Event_point_3 z_;
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include "Irrational_cross_section_impl.h"
#endif
#endif
