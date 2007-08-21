#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE



CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Irrational_cross_section CGAL_AOS3_TARG ::CS::Face_handle 
Irrational_cross_section CGAL_AOS3_TARG ::locate( CGAL_AOS3_TYPENAME Traits::Sphere_3_key  k) {

  /*if (ep.compare(z_, sweep_coordinate()) != CGAL::EQUAL) {
    std::cerr << z_ << std::endl;
    std::cerr << ep << std::endl;
    }
    CGAL_precondition(ep.compare(z_, sweep_coordinate()) == CGAL::EQUAL);*/
  std::vector<CGAL_AOS3_TYPENAME CS::Face_handle> faces;
  for (CGAL_AOS3_TYPENAME CS::Face_iterator it = cs_.faces_begin();it != cs_.faces_end(); ++it){
    faces.push_back(it);
  }
  return locate(faces.begin(), faces.end(), k);
}


CGAL_AOS3_TEMPLATE
bool 
Irrational_cross_section CGAL_AOS3_TARG ::locate_point_check_face(const CGAL_AOS3_TYPENAME Traits::Event_point_3& k,
								  CGAL_AOS3_TYPENAME CS::Face_const_handle it,
								  std::vector<int> &locations) const {
  CGAL_AOS3_TYPENAME CS::Halfedge_const_handle h= it->halfedge();
  //bool finite=false;
  do {
    //if (h->curve().is_finite()) {
    // we are outside
    if (!h->curve().is_finite()) {
      if (h->curve().is_right() && h->curve().is_inside()) return false;
    } else {
      CGAL_AOS3_TYPENAME Traits::Sphere_3_key sphere= h->curve().key();
      point_sphere_orientation(k, sphere, locations);
      //DPRINT(std::cout << "Testing " << h->curve() <<std::endl);
      if (!h->curve().is_compatible_location(locations[sphere.input_index()])) {
	//DPRINT(std::cout << "Nixed by edge " << h->curve() << std::endl);
	return false;
      } 
    }
    //finite=true;
    //} else {
    //std::cout << "Skipping infinite " << h->curve() << std::endl;
    //}
    h= h->next();
  } while (h != it->halfedge());
  /*if (finite) {
    std::cout << "Face " << std::endl;;
    do {
    std::cout << h->curve() << std::endl;
    h= h->next();
    } while (h != it->halfedge());
    std::cout << std::endl;
    }*/

  return true;
}





// cache checked arcs



CGAL_AOS3_TEMPLATE
bool 
Irrational_cross_section CGAL_AOS3_TARG ::locate_point_check_face_vertices(const CGAL_AOS3_TYPENAME Traits::Event_point_3 &k,
									   CGAL_AOS3_TYPENAME CS::Face_const_handle it) const {
  
  CGAL_AOS3_TYPENAME CS::Halfedge_const_handle h= it->halfedge();
  do {
    if (h->vertex()->point().is_sphere_sphere() 
	&& h->curve().key() != h->next()->curve().key()
	&& (!h->curve().is_inside() && !h->next()->curve().is_inside())) {
      // NOTE what about degeneracies?  not sure if I need to handle
      // them here
      CGAL_AOS3_TYPENAME CS::Point npt= h->vertex()->point();
      if (tr_.oriented_side_of_separating_plane(npt.sphere_key(0),
						npt.sphere_key(1),
						k) == CGAL::ON_NEGATIVE_SIDE) {
	//CGAL_assertion(0);
	std::cout << "Face nixed by vertex " << npt << std::endl;
	return false;
      }
    }
    h= h->next();
  } while (h != it->halfedge());
  
  return true;
}


CGAL_AOS3_TEMPLATE
int 
Irrational_cross_section CGAL_AOS3_TARG ::sphere_location(CGAL_AOS3_TYPENAME Traits::Sphere_point_3 p,
							  CGAL_AOS3_TYPENAME Traits::Sphere_3_key s) const {
  int r=0;
  CGAL::Bounded_side bs=tr_.bounded_side_of_sphere(s,
						   p);
    
  /*T::Event_point_3 ep=sphere_start(locate_point);
    T::Line_3 l= ep.line();
    T::Event_point_3 a(sphere(s), l);
    if (a.is_valid()) {
    T::Event_point_3 b(sphere(s), l.opposite());
    CGAL::Comparison_result ca= a.compare(ep, 2);
    CGAL::Comparison_result cb= b.compare(ep, 2);*/
    
  if (bs== CGAL::ON_BOUNDED_SIDE) {
    r |= CS::Curve::lIN_BIT;
  } else if (bs==CGAL::ON_BOUNDARY) {
    r |= CS::Curve::lIN_BIT;
    r |= CS::Curve::lOUT_BIT;
  } else {
    r |= CS::Curve::lOUT_BIT;
  }

  /*
    NOTE can make it more specific for intersection point location
  */
    
  CGAL::Comparison_result xo= tr_.compare_sphere_center_c(s, p,
							  plane_coordinate(0));
  CGAL::Comparison_result yo= tr_.compare_sphere_center_c(s,p, 
							  plane_coordinate(1));
  if (xo  != CGAL::LARGER) {
    r |= CS::Curve::lR_BIT;
  } 
  if (xo != CGAL::SMALLER){
    r |= CS::Curve::lL_BIT;
  }
  if (yo != CGAL::LARGER) {
    r |= CS::Curve::lT_BIT;
  } 
  if (yo != CGAL::SMALLER) {
    r |= CS::Curve::lB_BIT;
  }
  return r;
}




 



CGAL_AOS3_TEMPLATE
bool 
Irrational_cross_section CGAL_AOS3_TARG ::behind_arc(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
						     CGAL_AOS3_TYPENAME CS::Curve arc,
						     int location) const{
  CGAL_precondition(!arc.is_inside());
  CGAL_precondition(arc.is_compatible_location(location));
  CGAL_AOS3_TYPENAME Traits::Coordinate_index C=arc.is_weakly_incompatible(location);

  if (C.is_valid()) { // we know it is opposite the center
    /*CGAL::Comparison_result cr = t_.compare_sphere_center_c(arc.key(),
      ep, C);
      if (cr == CGAL::SMALLER && arc.is_negative()
      || cr == CGAL::LARGER && !arc.is_negative()) return false;*/
    bool b= tr_.start_point_in_slab_c(arc.key(), k, other_plane_coordinate(C));
    return b;
  } else {
    return false;
  }
}


CGAL_AOS3_TEMPLATE
bool 
Irrational_cross_section CGAL_AOS3_TARG ::behind_arc(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3& k,
						     CGAL_AOS3_TYPENAME CS::Curve arc,
						     int location) const{
  CGAL_precondition(!arc.is_inside());
  CGAL_precondition(arc.is_compatible_location(location));
  CGAL_AOS3_TYPENAME Traits::Coordinate_index C=arc.is_weakly_incompatible(location);

  if (C.is_valid()) { // we know it is opposite the center
    /*CGAL::Comparison_result cr = t_.compare_sphere_center_c(arc.key(),
      ep, C);
      if (cr == CGAL::SMALLER && arc.is_negative()
      || cr == CGAL::LARGER && !arc.is_negative()) return false;*/
    bool b= tr_.point_in_slab_c(arc.key(), k, other_plane_coordinate(C));
    return b;
  } else {
    return false;
  }
}

CGAL_AOS3_TEMPLATE
void 
Irrational_cross_section CGAL_AOS3_TARG ::point_sphere_orientation(CGAL_AOS3_TYPENAME Traits::Sphere_point_3 k,
								   CGAL_AOS3_TYPENAME Traits::Sphere_3_key sphere,
								   std::vector<int> &locations
								   /*,
								     std::vector<CS::Curve> &edges*/) const {
  if (locations[sphere.input_index()]==0){
    locations[sphere.input_index()]=sphere_location(k, sphere);
    //DPRINT(std::cout << "For sphere " << sphere << " got " << SL::decode(locations[sphere]) << std::endl);
  }
}
  






CGAL_AOS3_TEMPLATE
bool 
Irrational_cross_section CGAL_AOS3_TARG ::rule_shoot_compare_if_rational(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & , 
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










CGAL_AOS3_TEMPLATE
void 
Irrational_cross_section CGAL_AOS3_TARG ::debug_rule_shoot_check(CGAL::Comparison_result check, 
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





CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Irrational_cross_section CGAL_AOS3_TARG ::debug_rule_shoot_answer(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t, 
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


CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Irrational_cross_section CGAL_AOS3_TARG ::rule_shoot_compare_SR(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t, 
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


CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Irrational_cross_section CGAL_AOS3_TARG ::rule_shoot_compare_SS(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t, 
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




// return comparison of point on edge of face to the shot rule on the C coordinate
// i.e. SMALLER if the point on edge is smaller than the rule coordinate
CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Irrational_cross_section CGAL_AOS3_TARG ::rule_shoot_edge_vertex(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t,
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





CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Irrational_cross_section CGAL_AOS3_TARG ::CS::Halfedge_handle 
Irrational_cross_section CGAL_AOS3_TARG ::shoot_rule(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 & t,
						     CGAL_AOS3_TYPENAME CS::Face_handle f,
						     const CGAL_AOS3_TYPENAME Traits::Sphere_3_key & pt,
						     CGAL_AOS3_TYPENAME CS::Rule_direction rd) {
  //CGAL_assertion(t.compare(pt, CGAL_AOS3_INTERNAL_NS::Sweep_coordinate::object())== CGAL::EQUAL);
  // CGAL_assertion(pt.compare(z_, CGAL_AOS3_INTERNAL_NS::Sweep_coordinate::object())== CGAL::EQUAL);
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
							 tr_.sphere_events(pt).first, rd,
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





CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Irrational_cross_section CGAL_AOS3_TARG ::CS::Halfedge_handle
Irrational_cross_section CGAL_AOS3_TARG ::extend_rule(const CGAL_AOS3_TYPENAME Traits::Event_point_3 &t,
	    CGAL_AOS3_TYPENAME CS::Halfedge_handle h,
	    CGAL_AOS3_TYPENAME CS::Curve rule) {
  CGAL_AOS3_TYPENAME CS::Halfedge_handle v= find_rule_vertex(t,h->face(),rule);
  return cs_.split_face(rule, h, v);
}

CGAL_AOS3_TEMPLATE
void 
Irrational_cross_section CGAL_AOS3_TARG ::roll_back_rule(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &t,
	       CGAL_AOS3_TYPENAME CS::Halfedge_handle cur) {
  std::cout << "Rolling back, rolling back..." << cur->curve()
	    << " at " << cur->vertex()->point() << std::endl;
  std::vector<CGAL_AOS3_TYPENAME CS::Halfedge_handle> bits;
  bits.push_back(cur);
  do {
    bits.push_back(cs_.next_edge_on_rule(bits.back()));
  } while (bits.back() != CGAL_AOS3_TYPENAME CS::Halfedge_handle());
  bits.pop_back();
  
  if (bits.size() > 1) {
    
    // each internal intersection is either 3 way with only one side or 4 way
    for (unsigned int i=0; i< bits.size()-1; ++i) {
      if (cs_.degree(bits[i]->vertex())==4) {
	// nothing to do
      } else {
	std::cout << "Filling in vertex from edge ";
	cs_.write(bits[i], std::cout) << " with ";
	CGAL_AOS3_TYPENAME CS::Halfedge_handle ce= cs_.cross_edge(bits[i]); 
	cs_.write(ce, std::cout) << std::endl;
	CGAL_AOS3_TYPENAME CS::Face_handle f=ce->prev()->opposite()->face();
	
	CGAL_AOS3_TYPENAME CS::Halfedge_handle rt= extend_rule(t, 
							       cs_.find_halfedge(bits[i]->vertex(), f),
							       ce->opposite()->curve());
	
	
	//std::cout << "Adding rule from ";
	//cs_.write(find_halfedge(bits[i]->vertex(),f), std::cout) << " to ";
	//cs_.write(rt) << std::endl;
	//cs_.split_face(ce->opposite_curve(), find_halfedge(bits[i]->vertex(),f),
	//      rt);				    
      }
    }
    if (bits.back()->vertex()->point().is_sphere_extremum()
	&& !bits.back()->next()->curve().is_inside()) {
      // need the back
    } else {
      CGAL_AOS3_TYPENAME CS::Vertex_handle vh= bits.back()->vertex();
      cs_.join_face(bits.back(), true);
      //clean_up_vertex(t, vh);
    }
    for (unsigned int i=1; i< bits.size()-1; ++i) {
      cs_.join_face(bits[i], true);
    }
  }
  
}							   




CGAL_AOS3_TEMPLATE
void 
Irrational_cross_section CGAL_AOS3_TARG ::clean_up_vertex(const CGAL_AOS3_TYPENAME Traits::Event_point_3 &t,
		     CGAL_AOS3_TYPENAME CS::Vertex_handle vh) {
  if (cs_.is_redundant(vh)) {
    cs_.remove_vertex(vh);
  }
  if (vh->point().is_sphere_extremum()) {
    CGAL_AOS3_TYPENAME Traits::Sphere_3_key k= vh->point().sphere_key();
    CGAL_AOS3_TYPENAME CS::Halfedge_handle rule, arc, c= vh->halfedge();
    do {
      if (c->curve().is_rule() && c->curve().key() == k) {
	rule= c;
      } 
      if (!c->curve().is_rule()) {
	arc= c;
      }
      c= c->next()->opposite();
    } while (c != vh->halfedge());
    if (rule == CGAL_AOS3_TYPENAME CS::Halfedge_handle()) {
      if (arc->curve().is_inside()) arc=arc->opposite();
      CGAL_AOS3_TYPENAME CS::Curve rule= CS::Curve::make_rule(arc->curve().key(),
							      vh->point().sphere_extremum_index());
      CGAL_AOS3_TYPENAME CS::Halfedge_handle ovh
	= find_rule_vertex(t, arc->face(),
			   rule);
      cs_.split_face(rule, cs_.find_halfedge(vh, arc->face()), ovh);
    }
  }
  
}


CGAL_AOS3_TEMPLATE
 CGAL_AOS3_TYPENAME Irrational_cross_section CGAL_AOS3_TARG ::CS::Halfedge_handle 
 Irrational_cross_section CGAL_AOS3_TARG :: find_rule_vertex(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &t, 
		   CGAL_AOS3_TYPENAME CS::Face_handle f,
		   CGAL_AOS3_TYPENAME CS::Curve rule) {
    CGAL_AOS3_TYPENAME CS::Halfedge_handle h;
    std::cout << "Searching for vertex for rule " << rule << " in face ";
    cs_.write(f, std::cout) << std::endl;
    try {
      h= shoot_rule(t, f, rule.key(), rule.rule_direction());
      
      //check_edge_collapse(h->prev());
    } catch (On_vertex_exception e) {
      CGAL_AOS3_TYPENAME CS::Vertex_handle v= e.vertex_handle();
      // if it is a rule in the same direction, return it, otherwise
      // pick a random edge


      if (v->point().is_rule_rule() 
	  || v->point().is_sphere_rule() 
	  && v->point().rule_coordinate()== rule.constant_coordinate()) {
	// insert on vertex
	return cs_.find_halfedge(v,f);
      } else {
	// if I am shooting up, make sure I am above the point etc.
	/*Halfedge_handle h0= sds_.find_halfedge(v,f);
	  Halfedge_handle h1= h->next();
	  if (h0->curve() == h1->curve()) {
	  bool cum=false;
	  if (rule.is_vertical()) cum = !cum;
	  if (h0->curve().arc_index() ==0 || h0->curve().arc_index() ==2) cum= !cum;
	  if (cum) {
	  h= h1;
	  } else {
	  h= h0;
	  }
	  
	  } else {
	  CGAL_assertion(0);
	  }*/
	CGAL_AOS3_TYPENAME CS::Halfedge_handle vh= v->halfedge();
	do {
	  if ( vh->curve().is_arc() && vh->face() == f) {
	    h=vh;
	    break;
	  } else if (vh->curve().is_arc() && vh->opposite()->face() == f) {
	    h= vh->opposite();
	    break;
	  }
	  vh= vh->next()->opposite();
	} while (true);
      }
    }
    //cs_.clean_edge(h);
    CGAL_AOS3_TYPENAME CS::Vertex_handle nvh= cs_.new_vertex(CGAL_AOS3_TYPENAME CS::Point(rule, 
											  h->curve()));
    CGAL_AOS3_TYPENAME CS::Halfedge_handle v= cs_.insert_vertex(nvh, h);
    //check_edge_collapse(h);
    //CGAL_assertion(h->prev()->vertex()==v);
    CGAL_assertion(v->face() == f);
    return v;
  }


CGAL_AOS3_END_INTERNAL_NAMESPACE
