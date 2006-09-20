#include <CGAL/Arrangement_of_spheres_3/Slice.h>

#define DPRINT(x) x



Slice::Halfedge_handle Slice::shoot_rule(const T::Sphere_point_3 & t,
					 Face_handle f,
					 const T::Sphere_point_3 & pt,
					 Rule_direction rd) {
  //Sds::Curve rule= Sds::Curve::make_rule(ruledir);
  bool backwards= rd.is_backwards();
  DPRINT(std::cout << "Rule is " << rd << std::endl);
  
  Halfedge_handle end= f->halfedge();
  
  while (!rd.can_intersect(end->curve())) {
    end=end->next();
  }
  while (rd.can_intersect(end->curve())) {
    end=end->next();
  }
  Halfedge_handle h=end;
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
  return Halfedge_handle();
}



// return comparison of point on edge of face to the shot rule on the C coordinate
// i.e. SMALLER if the point on edge is smaller than the rule coordinate
CGAL::Comparison_result Slice::rule_shoot_edge_vertex(const T::Sphere_point_3 & t,
						      const T::Sphere_point_3 & pt, 
						      Rule_direction rd,
						      Sds::Curve hp,
						      Sds::Point p,
						      Sds::Curve hn) const {
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
				   (hn.quadrant() & Sds::Curve::T_BIT) 
				   || rd.is_vertical() && 
				   (hn.quadrant() & Sds::Curve::R_BIT));
    } else if (hp.is_arc() && hn.is_rule()){
      return rule_shoot_compare_SR(t, pt, rd, hp, hn.key(), 
				   p,
				   !rd.is_vertical() && 
				   (hp.quadrant() & Sds::Curve::T_BIT) 
				   || rd.is_vertical() && 
				   (hp.quadrant() & Sds::Curve::R_BIT));
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
      std::pair<T::Event_point_3, T::Event_point_3> ep3
	= t_.intersection_2_events(p.sphere_key(0),
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


CGAL::Comparison_result Slice::rule_shoot_compare_SS(const T::Sphere_point_3 & t, 
						     const T::Sphere_point_3 & pt, 
						     Rule_direction rd,
						     Sds::Point point) const {
  bool exact;
  CGAL::Comparison_result debug_answer= debug_rule_shoot_answer(t, pt, rd,
								point,
								exact);
  //std::cout << "SS " << arc0 << "--" << pt << "--" << arc1 << std::endl;
   
  T::Coordinate_index C= rd.constant_coordinate();

  CGAL::Comparison_result cr= t_.compare_sphere_sphere_at_sweep(t,
								point.sphere_key(0),
								point.sphere_key(1),
								pt, C);
  debug_rule_shoot_check(debug_answer, cr, exact);
  return cr;
}



CGAL::Comparison_result Slice::rule_shoot_compare_SR(const T::Sphere_point_3 & t, 
						     const T::Sphere_point_3 & pt,
						     Rule_direction rd,
						     Sds::Curve arc,
						     T::Key orule,
						     Sds::Point debug_pt,
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

  CGAL::Comparison_result c= t_.compare_sphere_center_c(arc.key(), pt,
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
    intersection_side=t_.bounded_side_of_sphere_projected(t, arc.key(),
							  orule,
							  pt,
							  rd.constant_coordinate());
  } else {
    intersection_side=t_.bounded_side_of_sphere_projected(t, arc.key(),
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


CGAL::Comparison_result Slice::debug_rule_shoot_answer(const T::Sphere_point_3 & t, 
						       const T::Sphere_point_3 & pt,
						       Rule_direction rd,
						       Sds::Point point,
						       bool &exact) const {
  //CGAL::Comparison_result answer;
  NT z;
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

  //Sds::Point pt(p,n);
  const T::Sphere_point_3 & sp= sphere_point_rz(point, z);
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


bool Slice::rule_shoot_compare_if_rational(const T::Sphere_point_3 & t, 
					   const T::Sphere_point_3 & pt,
					   Rule_direction rd,
					   Sds::Point point,
					   CGAL::Comparison_result &ret) const {
  //NT coord;
  if (point.is_rule_rule()) {
    ret= t_.compare_sphere_center_c(point.rule_key(rd.constant_coordinate()), pt,
				    rd.constant_coordinate());
    return true;
  } else if (point.is_sphere_rule() 
	     && point.rule_coordinate() == rd.constant_coordinate()) {
    ret= t_.compare_sphere_center_c(point.rule_key(), pt,
				    rd.constant_coordinate());
    return true;
  } else {
    return false;
  }
}

  

