#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section_rules.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE




#define DPRINT(x)







CGAL_AOS3_TEMPLATE
void 
Irrational_cross_section_rules CGAL_AOS3_TARG ::debug_check(CGAL::Comparison_result check, 
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
Irrational_cross_section_rules CGAL_AOS3_TARG ::debug_compare_to_answer(const Center_point_3 & pt,
									CGAL_AOS3_TYPENAME CS::Rule_direction rd,
									CGAL_AOS3_TYPENAME CS::Point point,
									bool &exact) const {
  //CGAL::Comparison_result answer;
  const Sphere_point_3 &t=pt.coord();
  Sphere_3_key center=pt.key();
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
  
  CGAL::Comparison_result cr=tr_.compare_to_rule_c(sp, center, rd.constant_coordinate());

  //sp.compare(pt, rd.constant_coordinate());
   
  DPRINT(std::cout << "Got " << cr << std::endl);
  if (cr== CGAL::EQUAL) {
    /*std::cout << "EQUAL for " << p << " " << n << ": " << sp 
      <<  " vs " << ep << " on " 
      << rule.constant_coordinate() << std::endl;*/
  }
  return Comparison_result(-cr);;
}










CGAL_AOS3_TEMPLATE
CGAL::Comparison_result 
Irrational_cross_section_rules CGAL_AOS3_TARG ::debug_compare_to_answer(const Sphere_point_3 & t,
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
  
  CGAL::Comparison_result cr=tr_.compare_c(sp, t, rd.constant_coordinate());

  //sp.compare(pt, rd.constant_coordinate());
   
  DPRINT(std::cout << "Got " << cr << std::endl);
  if (cr== CGAL::EQUAL) {
    /*std::cout << "EQUAL for " << p << " " << n << ": " << sp 
      <<  " vs " << ep << " on " 
      << rule.constant_coordinate() << std::endl;*/
  }
  return Comparison_result(-cr);;
}






















CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Irrational_cross_section_rules CGAL_AOS3_TARG ::CS::Halfedge_handle
Irrational_cross_section_rules CGAL_AOS3_TARG ::extend_rule(const CGAL_AOS3_TYPENAME Traits::Event_point_3 &t,
							    CGAL_AOS3_TYPENAME CS::Halfedge_handle h,
							    CGAL_AOS3_TYPENAME CS::Curve rule) {
  CGAL_AOS3_TYPENAME CS::Halfedge_handle v= find_rule_vertex(t,h->face(),rule);
  return cs_.split_face(rule, h, v);
}










CGAL_AOS3_TEMPLATE
void 
Irrational_cross_section_rules CGAL_AOS3_TARG ::roll_back_rule(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &t,
							       CGAL_AOS3_TYPENAME CS::Halfedge_handle cur) {
  CGAL_LOG(Log::LOTS, "Rolling back, rolling back..." << cur->curve()
	   << " at " << cur->vertex()->point() << std::endl);
  std::vector<CGAL_AOS3_TYPENAME CS::Halfedge_handle> bits;
  bits.push_back(cur);
  do {
    bits.push_back(cs_.next_edge_on_rule(bits.back()));
  } while (bits.back() != CGAL_AOS3_TYPENAME CS::Halfedge_handle());
  bits.pop_back();
  for (unsigned int i=0; i< bits.size(); ++i) {
    CGAL_LOG_WRITE(Log::LOTS, cs_.write(bits[i], LOG_STREAM) << ", ");
  }
  CGAL_LOG(Log::LOTS, std::endl);
  
  if (bits.size() > 1) {
    
    // each internal intersection is either 3 way with only one side or 4 way
    for (unsigned int i=0; i< bits.size()-1; ++i) {
      if (cs_.degree(bits[i]->vertex())==4) {
	// nothing to do
      } else {
	CGAL_LOG(Log::LOTS,"Filling in vertex from edge ");
	CGAL_LOG_WRITE(Log::LOTS, cs_.write(bits[i], LOG_STREAM) << " with ");
	CGAL_AOS3_TYPENAME CS::Halfedge_handle ce= cs_.cross_edge(bits[i]); 
	CGAL_LOG_WRITE(Log::LOTS, cs_.write(ce, LOG_STREAM) << std::endl);
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
Irrational_cross_section_rules CGAL_AOS3_TARG ::clean_up_vertex(const CGAL_AOS3_TYPENAME Traits::Event_point_3 &t,
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
CGAL_AOS3_TYPENAME Irrational_cross_section_rules CGAL_AOS3_TARG ::CS::Halfedge_handle 
Irrational_cross_section_rules CGAL_AOS3_TARG :: find_rule_vertex(const CGAL_AOS3_TYPENAME Traits::Sphere_point_3 &t, 
								  CGAL_AOS3_TYPENAME CS::Face_handle f,
								  CGAL_AOS3_TYPENAME CS::Curve rule) {
  CGAL_AOS3_TYPENAME CS::Halfedge_handle h;
  CGAL_LOG(Log::LOTS, "Searching for vertex for rule " << rule << " in face ");
  CGAL_LOG_WRITE(Log::LOTS, cs_.write(f, LOG_STREAM) << std::endl);
  try {
    h= shoot_rule(CGAL_AOS3_TYPENAME Traits::Center_point_3(rule.key(), t), f,  rule.rule_outward_direction());
      
    //check_edge_collapse(h->prev());
  } catch (On_vertex_exception e) {
    CGAL_AOS3_TYPENAME CS::Vertex_handle v= e.vertex_handle();
    // if it is a rule in the same direction, return it, otherwise
    // pick a random edge

    CGAL_LOG(Log::LOTS,  "rule shooting landed on vertex ");
    CGAL_LOG_WRITE(Log::LOTS, cs_.write(v, LOG_STREAM) << std::endl);

    if (v->point().is_rule_rule() 
	|| v->point().is_sphere_rule() 
	&& v->point().rule_constant_coordinate()== rule.constant_coordinate()) {
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
	CGAL_LOG(Log::LOTS,  "Inspecting ");
	CGAL_LOG_WRITE(Log::LOTS, cs_.write(vh, LOG_STREAM) << std::endl);
	if ( vh->curve().is_arc() && vh->face() == f && vh->curve().key() != rule.key()) {
	  h=vh;
	  break;
	} else if (vh->curve().is_arc() && vh->opposite()->face() == f && vh->curve().key() != rule.key()) {
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
