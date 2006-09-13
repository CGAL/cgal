#include <CGAL/Arrangement_of_spheres_3/Slice.h>

#define DPRINT(x)


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
    
  Halfedge_handle h= f->halfedge();
  std::vector<Halfedge_handle> edges;

  bool has_last=false;
  CGAL::Comparison_result last;
  if (rule.can_intersect(h->curve()) 
      && rule.can_intersect(h->prev()->curve())) {
    has_last=true;
    last= rule_shoot_edge_vertex(source, rule, h->prev()->curve(),
				 h->prev()->vertex()->point(),
				 h->curve());
  } else {
    DPRINT(std::cout << "Eliminating " <<  h->prev()->vertex()->point() << std::endl);
  }
  do {
    // later skip second one if first is wrong
    CGAL::Comparison_result pbs, nbs;
    DPRINT(std::cout << "\nTesting edge: " << h->prev()->vertex()->point() << "--" 
	      << h->curve() << "--"
	   << h->vertex()->point() << std::endl;)
    if (rule.can_intersect(h->curve())) {
      if (has_last) {
	DPRINT(std::cout << "Using last of " << last << std::endl);
	pbs=last;
      } else if (backwards) {
	pbs= CGAL::LARGER;
      } else {
	pbs=CGAL::SMALLER;
      }
      /* pbs= shoot_rule_edge_vertex(ep, rule, h->prev()->curve(),
	 h->prev()->vertex()->point(),
	 h->curve());
	 }*/
      if (rule.can_intersect(h->next()->curve())) {
	nbs= rule_shoot_edge_vertex(source, rule, h->curve(),
				    h->vertex()->point(),
				    h->next()->curve());
	DPRINT(std::cout << "Results are " << nbs << " " << pbs << std::endl);
	has_last=true;
	last = nbs;
      } else {
	has_last=false;
	if (backwards) nbs = CGAL::SMALLER;
	else nbs = CGAL::LARGER;
      }

      if (!backwards && pbs != CGAL::LARGER && nbs != CGAL::SMALLER
	  || backwards && pbs != CGAL::SMALLER && nbs != CGAL::LARGER) {
	//std::cout << "it passes" << std::endl;
	edges.push_back(h);
      }
    } else {
      //std::cout << "Filtered" <<  std::endl;
      has_last=false;
    }

      
    h= h->next();
    CGAL_assertion(h->face()== f);
  } while (h != f->halfedge());


  if (edges.size() == 0) {
    std::cerr << "No edge found " << std::endl;
    return Halfedge_handle();
  } else if (edges.size() ==1) {
    return edges[0];
  } else {

    CGAL_assertion(edges.size()== 2);
    std::cout << "Both edges " << edges[0]->curve() << " and " 
	      << edges[1]->curve()
	      << " are OK" << std::endl;
    if (edges[0]->vertex()== edges[1]->opposite()->vertex()) {
      throw On_vertex_exception(edges[0]->vertex());
    } else {
      CGAL_assertion(edges[1]->vertex()== edges[0]->opposite()->vertex());
      throw On_vertex_exception(edges[1]->vertex());
    }
  }
}



// return comparison of separator to intersection point on the C coordinate
// i.e. SMALLER if the separator is SMALLER than the intersection point
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
	// rule from same sphere
	CGAL_assertion(hp.is_inside());
	bool exact;
	CGAL::Comparison_result debug_answer= debug_rule_shoot_answer(ep, rule, 
								      Sds::Point(hp, hn), 
								      exact);
	if (hp.part() & hn.part() & rule.part()) {
	  // check rational should get this
	  CGAL_assertion(0);
	  /*CGAL::Comparison_result ret
	    = CGAL::compare(t_.center_c(hn.key(), rule.constant_coordinate()), 
			    ep.simple_coordinate(rule.constant_coordinate()));
			    debug_rule_shoot_check(debug_answer, ret, exact);*/
	  return CGAL::SMALLER;
	}  else if (hp.is_top()&& hn.is_top()
		    || hp.is_right() && hn.is_right()
		    /*hp.part() & hn.part() & Sds::Curve::T_BIT
		      || hp.part() & hn.part() & Sds::Curve::R_BIT*/) {
	  debug_rule_shoot_check(debug_answer, CGAL::LARGER, exact);
	  return CGAL::LARGER;
	} else {
	  debug_rule_shoot_check(debug_answer, CGAL::SMALLER, exact);
	  return CGAL::SMALLER;
	}
      }
    }
  } else {
    return rule_shoot_compare_SS(ep, rule, hp, p, hn);
  }

  CGAL_assertion(0);
  return CGAL::SMALLER;
}


CGAL::Comparison_result Slice::rule_shoot_compare_SS(const T::Sphere_point_3 & ep, 
						     Sds::Curve srule,
						     Sds::Curve arc0,
						     Sds::Point pt,
						     Sds::Curve arc1) const {
  bool exact;
  CGAL::Comparison_result debug_answer= debug_rule_shoot_answer(ep, srule,
								pt,
								exact);
  //std::cout << "SS " << arc0 << "--" << pt << "--" << arc1 << std::endl;
   
  T::Coordinate_index C= srule.constant_coordinate();
  bool tangent= (arc0.key()==arc1.key());
  if (tangent) {
    CGAL_assertion(0);
    std::cout << "Tangent." << std::endl;
    /*if (arc0.key() == pt.sphere_key(0)) arc1= pt.sphere(1);
    else arc1= pt.sphere(0);
    std::cout << "Using " << arc0 << "--" << pt << "--" << arc1 << std::endl;*/
  }
  // NOTE maybe make these into index based predicates
  bool hi= t_.sphere_intersects_rule(arc0.key(), srule.key(), C);
  bool ohi= t_.sphere_intersects_rule(arc0.key(), srule.key(), C);
  CGAL::Comparison_result hp= t_.compare_sphere_centers_c(arc0.key(),
							  srule.key(), C);
  CGAL::Comparison_result ohp= t_.compare_sphere_centers_c(arc1.key(), 
							   srule.key(), C);
  //std::cout << "Comparisons are " << hi << " " << hp << " and " << ohi << " " << ohp << std::endl;

  if (hi && !ohi || !ohi && !hi) {
    debug_rule_shoot_check(debug_answer,ohp, exact);
    return ohp;
  } else if (ohi && !hi) {
    debug_rule_shoot_check(debug_answer,hp, exact);
    return hp;
  } /*else if (!ohp && !hp) {
    debug_rule_shoot_check(debug_answer,hp, exact);
    CGAL_assertion(hp == ohp);
    return ohp;
    }*/
  // now they both intersect the plane


  // NOTE  hmmm, need to clean
  CGAL::Bounded_side bs
    = t_.bounded_side_of_sphere_on_equipower_plane_rule_line(arc0.key(), 
							      arc1.key(),
							      srule.key(),
							      T::Coordinate_index(C),
							     ep);
  if (bs == CGAL::ON_UNBOUNDED_SIDE){
    //std::cout << "Missed." << std::endl;
    CGAL::Comparison_result ret= t_.compare_equipower_point_to_rule(arc0.key(),
								 arc1.key(),
								 srule.key(), C);
    debug_rule_shoot_check(debug_answer,ret, exact);
    return ret;
  } else if (bs== CGAL::ON_BOUNDARY) {
    //std::cout << "Equal one of them " << std::endl;
    debug_rule_shoot_check(debug_answer,CGAL::EQUAL, exact);
    return CGAL::EQUAL;
  } else {
    if (arc0.is_inside() || arc1.is_inside()) {
      // arcs are at extremum of face if I am checking
      CGAL::Sign sn= t_.sign_of_separating_plane_normal_c(pt.sphere_key(0),
						       pt.sphere_key(1), C);
      CGAL_assertion(sn != CGAL::ZERO);
      CGAL::Comparison_result cret= CGAL::enum_cast<CGAL::Comparison_result>(sn);
      debug_rule_shoot_check(debug_answer, 
			     cret, 
			     exact);
      return cret;
    } else {
      //bool between =true;
      CGAL::Sign dir= t_.sign_of_equipower_plane_normal_c(arc0.key(), arc1.key(),
						       C);
      //CGAL::Comparison_result dir=CGAL::compare(eqp.orthogonal_vector()[C], 0); 
      int cum=0;
      CGAL::Oriented_side os= t_.oriented_side_of_equipower_plane(arc0.key(), arc1.key(), ep);
      if (os == CGAL::ON_POSITIVE_SIDE) ++cum;
      if (dir == CGAL::POSITIVE) ++cum;
      if (cum%2==0) {
	debug_rule_shoot_check(debug_answer,CGAL::SMALLER, exact);
	return CGAL::SMALLER;
      } else {
	debug_rule_shoot_check(debug_answer,CGAL::LARGER, exact);
	return CGAL::LARGER;
      }
    }
  }
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

   
      
  CGAL::Bounded_side intersection_side;
  if (srule.constant_coordinate() == plane_coordinate(0)) {
    intersection_side= t_.bounded_side_of_sphere(arc.key(),
						 srule.key(),
						 orule,
						 ep);
  } else {
    intersection_side= t_.bounded_side_of_sphere(arc.key(),
						 orule,
						 srule.key(),
						 ep);
  }
  //ep, orule, 
  //							       arc.key(),
  //							       srule.constant_coordinate());
  

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
    
  CGAL::Comparison_result ans= CGAL::Comparison_result(-t_.compare_sphere_center_c(rule.key(), sp,
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


bool Slice::rule_shoot_compare_if_rational_arc(const T::Sphere_point_3 & ep,
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
}


bool Slice::rule_shoot_compare_if_rational(const T::Sphere_point_3 & ep,
					   Sds::Curve rule,
					   Sds::Point pt,
					   CGAL::Comparison_result &ret) const {
  //NT coord;
  if (pt.is_rule_rule()) {
    ret= t_.compare_sphere_centers_c(pt.rule_key(rule.constant_coordinate()), 
				     rule.key(),
				     rule.constant_coordinate());
    return true;
  } else if (pt.is_sphere_rule() 
	     && pt.rule_coordinate() == rule.constant_coordinate()) {
    ret= t_.compare_sphere_centers_c(pt.rule_key(),
				       rule.key(),
				       rule.constant_coordinate());
    return true;
  } else {
    // need the vertex_handle
    // actually, I don't think it was ever useful
    return false;
  }
  /*if (a.is_arc() && rule_shoot_compare_if_rational_arc(ep, rule, a, ret)){
      return true;
    }
    if (b.is_arc() && rule_shoot_compare_if_rational_arc(ep, rule, b, ret)){
      return true;
    }
    
    return false;
    }*/
  /*std::cout << "It is rational with coordinate " 
    << CGAL::to_double(coord) << std::endl;*/

  /*ret= CGAL::compare(coord, 
		     t_.center(rule.key())[rule.constant_coordinate()]);
		     return true;*/
}

  

// compare along C where the intersection point of the two rules
// (defined by Coord C of ep) is relative to the sphere of the arc
#if 0
CGAL::Bounded_side Slice::rule_shoot_source_side(const T::Sphere_point_3 & ep,
						 T::Key rule,
						 Sds::Curve orule,
						 T::Key arc_index,
						 int C) const {
  return sc_.bounded_side_of_sphere_rule_rule_line(epi,
						   orule.index(),
						   arc_index,
						   C,
						   ep);
  /*NT p[2];
  p[C]= ep.simple_coordinate(C);
  p[1-C]= center_c(orule.key(), orule.constant_coordinate());
    
  std::cout << "Line through " << CGAL::to_double(p[0]) 
	    << " " << CGAL::to_double(p[1]) << std::endl;
  T::Line_3 l(T::Point_3(p[0], p[1], 0), T::Vector_3(0,0,1));
    
  const T::Sphere_point_3 & p0=const T::Sphere_point_3 &(sphere(arc_index), l);
  CGAL::Bounded_side intersection_side;
  if (!p0.is_valid()) {
    intersection_side=CGAL::ON_UNBOUNDED_SIDE;
  } else if (p0 == ep) {
    std::cout << "Equal" << std::endl;
    intersection_side=  CGAL::ON_BOUNDARY;
  } else {
    const T::Sphere_point_3 & p1=const T::Sphere_point_3 &(sphere(arc_index), 
					 l.opposite());
    if (p1 == ep) {
      std::cout << "Equal 1" << std::endl;
      intersection_side=  CGAL::ON_BOUNDARY;
    }
    CGAL_assertion(p0 <= p1);
    if (p0 < ep && p1 > ep){
      intersection_side= CGAL::ON_BOUNDED_SIDE;
    } else {
      intersection_side= CGAL::ON_UNBOUNDED_SIDE;
    }
  }
  return intersection_side;*/
}
#endif
