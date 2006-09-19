#include <CGAL/Arrangement_of_spheres_3/Slice.h>

#define DPRINT(x) x


// hits in order T R B L
Slice::Halfedge_handle Slice::shoot_rule(const T::Sphere_point_3& source,
					       Face_handle f,
					       int type) {
  //Temp_point tp(spheres_, source.sphere());
  t_.set_temp_sphere(source.sphere());
  return shoot_rule(source, f, Sds::Curve(T::Key::target_key(),  Sds::Curve::Part(type)));

}



Slice::Halfedge_handle Slice::shoot_rule(const T::Sphere_point_3 & source,
					 Face_handle f,
					 Sds::Curve rule) {
  bool backwards= (rule.is_top() || rule.is_left());
  DPRINT(std::cout << "Rule is " << rule << std::endl);
  
  Halfedge_handle end= f->halfedge();
  
  while (!rule.can_intersect(end->curve())) {
    end=end->next();
  }
  while (rule.can_intersect(end->curve())) {
    end=end->next();
  }
  Halfedge_handle h=end;
  while (!rule.can_intersect(h->curve())) {
    h=h->next();
  }
  
  do {
    // later skip second one if first is wrong
    DPRINT(sds_.write(h, std::cout) << " is being tested" << std::endl);
    if (h->next() == end) {
      DPRINT(sds_.write(h, std::cout) << " is returned by default" << std::endl);
      return h;
    } else {
      CGAL::Comparison_result cr= rule_shoot_edge_vertex(source, rule, h->curve(),
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
  return Halfedge_handle();
}



// return comparison of point on edge of face to the shot rule on the C coordinate
// i.e. SMALLER if the point on edge is smaller than the rule coordinate
CGAL::Comparison_result Slice::rule_shoot_edge_vertex(const T::Sphere_point_3 & ep, 
						      Sds::Curve rule,
						      Sds::Curve hp,
						      Sds::Point p,
						      Sds::Curve hn) const {
  //CGAL_assertion(h != b);
  //CGAL_assertion(h->curve() != b->curve());
  CGAL::Comparison_result ret;
  // add filters on sphere centers when arcs are involved.
  if (rule_shoot_compare_if_rational(ep, rule, p, ret)){
     
    bool exact;
    CGAL::Comparison_result debug_answer= debug_rule_shoot_answer(ep, rule,
								  p, 
								  exact);
    debug_rule_shoot_check(debug_answer, ret, exact);
    return ret;
  } else if (p.is_sphere_rule()) {
    if (hp.is_rule() && hn.is_arc()){
      return rule_shoot_compare_SR(ep, rule, hn, hp.key(),
				   p,
				   !rule.is_vertical() && 
				   (hn.quadrant() & Sds::Curve::T_BIT) 
				   || rule.is_vertical() && 
				   (hn.quadrant() & Sds::Curve::R_BIT));
    } else if (hp.is_arc() && hn.is_rule()){
      return rule_shoot_compare_SR(ep, rule, hp, hn.key(), 
				   p,
				   !rule.is_vertical() && 
				   (hp.quadrant() & Sds::Curve::T_BIT) 
				   || rule.is_vertical() && 
				   (hp.quadrant() & Sds::Curve::R_BIT));
    } else if (hp.key() == hn.key()){
      CGAL_assertion(hp.is_arc() && hn.is_arc());
      if (hp == hn) {
	//CGAL_assertion(hp==hn);
	bool arc_top;
	if (rule.constant_coordinate() == plane_coordinate(0)) {
	  arc_top= hp.is_right();
	} else {
	  arc_top = hp.is_top();
	}
	if (!hp.is_inside()) arc_top= !arc_top;
	  
	return rule_shoot_compare_SR(ep, rule, hp, p.rule_key(),p, arc_top);
      } else {
	// compare rational should have picked this up
	CGAL_assertion(0);
      }
    }
  } else {
    bool tangent= (hp.key()== hn.key());
    if (tangent) {
      std::pair<T::Event_point_3, T::Event_point_3> ep3
	= t_.intersection_2_events(p.sphere_key(0),
				   p.sphere_key(1));
      std::cout << "Tangent." << std::endl;
      CGAL::Comparison_result cr= ep3.first.compare(ep, sweep_coordinate());
      //debug_rule_shoot_check(debug_answer, cr, exact);
      return cr;
    } else {
      return rule_shoot_compare_SS(ep, rule,p);
    }
  }

  CGAL_assertion(0);
  return CGAL::SMALLER;
}


CGAL::Comparison_result Slice::rule_shoot_compare_SS(const T::Sphere_point_3 & ep, 
						     Sds::Curve srule,
						     Sds::Point pt) const {
  bool exact;
  CGAL::Comparison_result debug_answer= debug_rule_shoot_answer(ep, srule,
								pt,
								exact);
  //std::cout << "SS " << arc0 << "--" << pt << "--" << arc1 << std::endl;
   
  T::Coordinate_index C= srule.constant_coordinate();

  CGAL::Comparison_result cr= t_.compare_sphere_sphere_at_sweep(pt.sphere_key(0),
								pt.sphere_key(1),
								ep, ep, C);
  debug_rule_shoot_check(debug_answer, cr, exact);
  return cr;
}



CGAL::Comparison_result Slice::rule_shoot_compare_SR(const T::Sphere_point_3 & ep, 
						     Sds::Curve srule,
						     Sds::Curve arc,
						     T::Key orule,
						     Sds::Point debug_pt,
						     bool arc_above) const {
  bool exact;
  CGAL::Comparison_result debug_answer= debug_rule_shoot_answer(ep, srule, 
								debug_pt,
								exact);
  DPRINT(std::cout << "SR " << arc << "-" << arc_above << "-" << orule 
	 << ": " << srule << std::endl);
  //CGAL_precondition(srule.is_vertical() != orule.is_vertical());
  CGAL_precondition(arc.is_arc());
  //CGAL_precondition(orule.is_rule());
  //CGAL_precondition(srule.constant_coordinate() != orule.constant_coordinate());

  CGAL::Comparison_result c= t_.compare_sphere_center_c(arc.key(), ep,
							srule.constant_coordinate());
  if (srule.is_vertical()) {
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

      
  CGAL::Bounded_side intersection_side
    =t_.bounded_side_of_sphere_projected(arc.key(),
					 ep,
					 orule,
					 srule.constant_coordinate());
  

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


CGAL::Comparison_result Slice::debug_rule_shoot_answer(const T::Sphere_point_3 & ep, 
						       Sds::Curve rule,
						       /*Sds::Curve p,
							 Sds::Curve n,*/
						       Sds::Point pt,
						       bool &exact) const {
  //CGAL::Comparison_result answer;
  NT z;
  if (ep.has_simple_coordinate(sweep_coordinate())) {
    z= ep.simple_coordinate(sweep_coordinate());
    exact=true;
  } else {
    z= ep.approximate_coordinate(sweep_coordinate());
    exact=false;
  }
  /*DPRINT(std::cout << "Comparing " << t_.center(rule.key()) << " on coord " 
	 << rule.constant_coordinate() << " which is " 
	 << t_.center(rule.key())[rule.constant_coordinate()] << std::endl);*/

  //std::cout << "Z is " << z << " and exact is " << exact << std::endl;

  //Sds::Point pt(p,n);
  const T::Sphere_point_3 & sp= sphere_point_rz(pt, z);
  /*if (!sp.is_valid()) {
    const T::Sphere_point_3 & sp2= sphere_point_rz(pt, z);
    }*/
  CGAL_assertion(sp.is_valid());
    
  DPRINT(std::cout << "The point " << pt << " is " 
	 << CGAL::to_double(sp.exact_coordinate(plane_coordinate(0))) << " " 
	 << CGAL::to_double(sp.exact_coordinate(plane_coordinate(1))) << std::endl);
    
  CGAL::Comparison_result ans
    = CGAL::Comparison_result(-t_.compare_sphere_center_c(rule.key(), sp,
							  rule.constant_coordinate()));
  DPRINT(std::cout << "Got " << ans << std::endl);
  if (ans== CGAL::EQUAL) {
    /*std::cout << "EQUAL for " << p << " " << n << ": " << sp 
	      <<  " vs " << ep << " on " 
	      << rule.constant_coordinate() << std::endl;*/
  }
  return ans;
}

  
  
void Slice::debug_rule_shoot_check(CGAL::Comparison_result check, 
				   CGAL::Comparison_result computed,
				   bool exact) {
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


bool Slice::rule_shoot_compare_if_rational(const T::Sphere_point_3 & ep,
					   Sds::Curve rule,
					   Sds::Point pt,
					   CGAL::Comparison_result &ret) const {
  //NT coord;
  if (pt.is_rule_rule()) {
    ret= t_.compare_sphere_center_c(pt.rule_key(rule.constant_coordinate()), 
				    ep,
				    rule.constant_coordinate());
    return true;
  } else if (pt.is_sphere_rule() 
	     && pt.rule_coordinate() == rule.constant_coordinate()) {
    ret= t_.compare_sphere_center_c(pt.rule_key(),
				    ep,
				    rule.constant_coordinate());
    return true;
  } else {
    return false;
  }
}

  

