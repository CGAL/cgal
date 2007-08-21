#include <CGAL/Arrangement_of_spheres_3/Event_processor.h>
#include <algorithm>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE



CGAL_AOS3_TEMPLATE
Event_processor CGAL_AOS3_TARG::Event_processor(CCS &cs, Traits &tr): cs_(cs), tr_(tr){
  
}



CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::insert(Sphere_3_key k) {  
  ICSI icsi(tr_, cs_);
  try {
    CGAL_AOS3_TYPENAME CCS::Face_handle f= icsi.locate(k);
    //slice.new_marked_face(f);
    icsi.insert(k, f);
  } catch (CGAL_AOS3_TYPENAME ICSI::On_edge_exception e) {
    icsi.insert(k, e.halfedge_handle());
  } catch (CGAL_AOS3_TYPENAME ICSI::On_vertex_exception v) {
    icsi.insert(k, v.vertex_handle());
  }
}

CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::remove(Sphere_3_key k) {
  ICSR icsr(tr_, cs_);
  icsr.remove_sphere(k);
}

CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::intersect(Sphere_3_key k, Sphere_3_key l) {
  std::vector<CGAL_AOS3_TYPENAME CCS::Face_handle> border_faces_k,
    border_faces_l, border_both;
  // can I determine if I need inside or outside or both
  // do point/sphere location for insertion and use that to filter further
  gather_incident_faces(k, std::back_inserter(border_faces_k));
  gather_incident_faces(l, std::back_inserter(border_faces_l));
  std::set_intersection(border_faces_k.begin(), border_faces_k.end(),
			border_faces_l.begin(), border_faces_l.end(),
			std::back_inserter(border_both),
			CGAL_AOS3_TYPENAME CCS::Handle_compare());
  CGAL_AOS3_TYPENAME CCS::Face_handle f;
  if (border_both.size() > 1) {
    ICS ics(tr_, cs_);
    
    f=ics.locate(border_both.begin(), border_both.end(),
		 cs_.visitor().simulator()->current_time());
    //f= ics.locate(border_both.begin(), border_both.end());
  } else {
    f= border_both.back();
  }

  std::cout << "Intersecting on face ";
  cs_.write(f, std::cout) << std::endl;

  Halfedge_handle hk=f->halfedge();
  do {
    hk=hk->next();
  } while (!hk->curve().is_arc() || hk->curve().key() != k);
  Halfedge_handle hl=f->halfedge();
  do {
    hl=hl->next();
  } while (!hl->curve().is_arc() || hl->curve().key() != l);

  std::cout << "Intersection edges are ";
  cs_.write(hk, std::cout) << " and ";
  cs_.write(hl, std::cout) << std::endl;

  cs_.intersect_arcs(hl, hk);
}

CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::unintersect(Sphere_3_key k, Sphere_3_key l) {
  CGAL_AOS3_TYPENAME CCS::Halfedge_handle hl, hk= cs_.a_halfedge(k)->opposite();
  CGAL_AOS3_TYPENAME CCS::Halfedge_handle khs=hk;
  do {
    if (hk->next()->curve().key() == l 
	&& hk->next()->next() == hk){
      hl= hk->next();
      break;
    } else if (hk->opposite()->next()->curve().key() ==l
	       && hk->opposite()->next()->next() == hk->opposite()) {
      hl= hk->opposite()->next();
      hk= hk->opposite();
      break;
    }
    hk= cs_.next_edge_on_arc(hk);
  } while (hk != khs);
  // find shared face, f, and edges hk and hl
  cs_.unintersect_arcs(hl, hk);
}


CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::intersect(Sphere_3_key k, Sphere_3_key l, Sphere_3_key m) {
  Halfedge_handle hl, hm, hk= cs_.a_halfedge(k)->opposite();
  Halfedge_handle khs=hk;
  do {
    if (hk->next()->curve().key() == l 
	&& hk->next()->next()->curve().key() == m
	&& hk->next()->next()->next() == hk){
      hl= hk->next();
      hm= hk->next()->next();
      break;
    } else if (hk->next()->curve().key() == m 
	       && hk->next()->next()->curve().key() == l
	       && hk->next()->next()->next() == hk){
      hl= hk->next()->next();
      hm= hk->next();
      break;
    } else if (hk->opposite()->next()->curve().key() == l 
	       && hk->opposite()->next()->next()->curve().key() == m
	&& hk->opposite()->next()->next()->next() == hk){
      hk=hk->opposite();
      hl= hk->next();
      hm= hk->next()->next();
      break;
    } else if (hk->opposite()->next()->curve().key() == m 
	       && hk->opposite()->next()->next()->curve().key() == l
	       && hk->opposite()->next()->next()->next() == hk){
      hk=hk->opposite();
      hl= hk->next()->next();
      hm= hk->next();
      break;
    } 
    hk= cs_.next_edge_on_arc(hk);
  } while (hk != khs);
 

  Halfedge_handle hn= cs_.next_edge_on_arc(hk)->opposite();
  Halfedge_handle hp= cs_.next_edge_on_arc(hk);
  
  Halfedge_handle hnt= hn->next()->next();
  Halfedge_handle hpt= hp->prev()->prev();

  Vertex_handle nvhn= cs_.new_vertex(Point(hk->curve(), hnt->curve()));
  Vertex_handle nvhp= cs_.new_vertex(Point(hk->curve(), hpt->curve()));

  cs_.insert_vertex(nvhn, hnt);
  cs_.insert_vertex(nvhp, hpt);
  cs_.move_target(hp, nvhp, true);
  cs_.move_target(hn, nvhn, true);


  cs_.split_face(hk->curve(), hn->next()->opposite(), 
		 hp->prev()->opposite()->prev());
  cs_.join_face(hk, true);
}

CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::unintersect(Sphere_3_key k, Sphere_3_key l, Sphere_3_key m) {
  intersect(k,l,m);
}


CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::process_aar(Halfedge_handle h) {
  clear_event(h);
  CGAL_precondition(h->vertex()->point().is_sphere_rule());
  Halfedge_handle rule= cs_.cross_edge(h);
  CGAL_assertion(rule->vertex() ==h->vertex());
  Halfedge_handle target;
  if (h->next() == rule->opposite()) {
    target= h->prev();
  } else {
    target=h->opposite()->next();
  }
  Point pt= h->vertex()->point();
  Vertex_handle nv= cs_.new_vertex(pt);
  cs_.insert_vertex(nv, target);
  cs_.move_target(rule, nv, true);
}



CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::process_rar(Halfedge_handle h) {
  clear_event(h);
  if (h->opposite()->vertex()->point().is_sphere_extremum()) {
    h=h->opposite();
  }
  CGAL_assertion(cs_.degree(h->vertex())==3);
  CGAL_assertion(cs_.degree(h->opposite()->vertex())==3);
  // pick which to move, make sure don't move extremum
  Halfedge_handle move_rule= cs_.cross_edge(h->opposite())->opposite();
  CGAL_assertion(move_rule->vertex() == h->opposite()->vertex());
  Halfedge_handle onto_edge= cs_.cross_edge(h);
  if (onto_edge->face() != move_rule->face() 
      && onto_edge->face() != move_rule->opposite()->face()){
    onto_edge= onto_edge->opposite();
  }
  Halfedge_handle other_rule= cs_.cross_edge(h);
  // handle them not being on the same side
  if (onto_edge->face() != move_rule->face() 
      && onto_edge->face() != move_rule->opposite()->face()){
    Vertex_handle ov= move_rule->vertex();
    cs_.move_target(move_rule, other_rule->vertex(), false);
    cs_.move_target(other_rule, ov, true);
  } else {
    Point npt(move_rule->curve(), other_rule->curve());
    Vertex_handle nv= cs_.new_vertex(npt);
    cs_.insert_vertex(nv, other_rule);
    cs_.move_target(move_rule, nv, true);
  }
}


CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::clear_event(Halfedge_handle h) {
  h->set_event(cs_.visitor().simulator()->null_event());
  h->opposite()->set_event(cs_.visitor().simulator()->null_event());
}
  

CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::process_arr(Halfedge_handle h) {
  clear_event(h);
  if (!h->vertex()->point().is_rule_rule()) h= h->opposite();
  CGAL_assertion(h->vertex()->point().is_rule_rule());
  Halfedge_handle xr= cs_.cross_edge(h)->opposite(); // inward pointing
  //CGAL_assertion(xr->vertex() == h->vertex());
  std::cout << "XR is ";
  cs_.write(xr, std::cout) << std::endl;

  ICS ics(tr_, cs_);
  ics.roll_back_rule(cs_.visitor().simulator()->current_time(), xr);
  Halfedge_handle nxr= cs_.next_edge_on_rule(xr);
  if (nxr != Halfedge_handle()) {
    if (cs_.is_redundant(xr)) {
      cs_.join_face(xr, true);
      xr= nxr->opposite();
      std::cout << "XR is ";
      cs_.write(xr, std::cout) << std::endl;
      nxr=Halfedge_handle();
    } else {
      nxr= nxr->opposite();
      std::cout << "NXR is ";
      cs_.write(nxr, std::cout) << std::endl;
    }
  } else {
    if (cs_.is_redundant(xr)) {
      cs_.join_face(xr, true);
    } else {
      do {
	Halfedge_handle target_h;
	if (h->face() == xr->face() || h->face() == xr->opposite()->face()) {
	  target_h= h->prev();
	} else {
	  target_h=h->opposite()->next();
	}
	std::cout << "target is ";
	cs_.write(target_h, std::cout) << std::endl;
	
	Point pt(xr->curve(), target_h->curve());
	Vertex_handle nv= cs_.new_vertex(pt);
	cs_.insert_vertex(nv, target_h);
	cs_.move_target(xr, nv, true);
      } while ((xr=nxr) != Halfedge_handle()); // do other side if needed
    }
  }
 
 
  // could fail if there is a degeneracy

 // now the cross rule should be a t pointing on direction;
  
 
}

CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::process_aae(Sphere_3_key k, Sphere_3_key l, int i) {
  int index=i+1;
  if (tr_.compare_sphere_center_c(k, cs_.visitor().simulator()->current_time(), plane_coordinate(1-i))== CGAL::LARGER) {
    index+=2;
  }
  std::cout << "Extremum event between " << k << " and " << l << " on rules ";
  if (i==0) std::cout << "TB";
  else std::cout << "LR";
  std::cout << "   " << index << std::endl;

  // find outwards pointing rule edge
  Halfedge_handle hk= cs_.a_halfedge(k);
  while (!hk->vertex()->point().is_sphere_extremum() 
	 || hk->vertex()->point().is_sphere_extremum()  
	 && hk->vertex()->point().sphere_extremum_index().index() != index) {
    cs_.write(hk, std::cout) << " " << hk->vertex()->point().is_sphere_extremum() 
			     << " " <<  (hk->vertex()->point().is_sphere_extremum()? hk->vertex()->point().sphere_extremum_index().index(): -1) << std::endl;
    hk= cs_.next_edge_on_arc(hk);
  }
  Halfedge_handle rule= hk->opposite()->prev()->opposite(); // outward rule
 

  std::cout << index << " hk is ";
  cs_.write(hk, std::cout) << " and rule is ";
  cs_.write(rule, std::cout) << std::endl;

  // check if it ends on arc of l
  if (rule->vertex()->point().is_sphere_rule() 
      && rule->vertex()->point().sphere_key() == l) {
    
    Halfedge_handle new_extremum; // inside arc
    Curve base_curve, rule_curve=rule->curve();
    Point pt= hk->vertex()->point();
    if (rule->opposite()->prev()->prev() != hk->opposite()) {
      new_extremum=cs_.next_edge_on_arc(cs_.next_edge_on_arc(hk));
      base_curve= hk->curve();
    } else {
      CGAL_assertion(rule->opposite()->prev()->prev() == hk->opposite());
      base_curve= hk->next()->curve();
      new_extremum=cs_.next_edge_on_arc(hk->opposite())->opposite();
    }
    
    std::cout << "New extremum is ";
    cs_.write(new_extremum, std::cout) << " and base curve is " << base_curve << std::endl;

    CGAL_postcondition(new_extremum->curve().is_inside());
    cs_.set_curve(hk, base_curve);
    cs_.set_curve(cs_.next_edge_on_arc(hk), base_curve);
    
    cs_.join_face(rule, true);
    //cs_.remove_vertex(v0, ->curve());
    //cs_.remove_vertex(v1);
    Vertex_handle nvh= cs_.new_vertex(pt);
    Halfedge_handle nh=cs_.insert_vertex(nvh, new_extremum);
    cs_.set_curve(nh->opposite(), cs_.next_edge_on_arc(nh->opposite())->curve());
    cs_.set_curve(nh->next(), cs_.next_edge_on_arc(nh->next())->curve());
    
    ICS ics(tr_, cs_);
    Face_handle f;
    if (nh->curve().is_inside()) {
      nh=nh->opposite()->prev();
    }
    f= nh->face();
    Halfedge_handle hd= ics.find_rule_vertex(cs_.visitor().simulator()->current_time(),
			 f,
			 rule_curve);
    cs_.split_face(rule_curve, nh, hd);
  } else {
    ICS ics(tr_, cs_);
    ics.roll_back_rule(cs_.visitor().simulator()->current_time(), rule);

    Halfedge_handle new_extremum, new_target; // inside arc
    Curve base_curve, rule_curve=rule->curve();
    Point pt= hk->vertex()->point();
    if (rule->prev()->prev()->curve().is_arc()
	&& rule->prev()->prev()->curve().key() == l) {
      new_extremum=cs_.next_edge_on_arc(cs_.next_edge_on_arc(hk));
      base_curve= hk->curve();
      new_target= rule->prev()->prev()->opposite();
    } else {
      base_curve= cs_.next_edge_on_arc(hk)->curve();
      new_extremum=cs_.next_edge_on_arc(hk->opposite())->opposite();
      new_target= rule->opposite()->next()->next()->opposite();
    }

    std::cout << "New extremum is ";
    cs_.write(new_extremum, std::cout) << " and base curve is " << base_curve 
			    << " and new target is ";
    cs_.write(new_target, std::cout) << std::endl;

    
    CGAL_assertion(new_target->curve().key() ==l);
    CGAL_postcondition(new_extremum->curve().is_inside());
    cs_.set_curve(hk, base_curve);
    cs_.set_curve(cs_.next_edge_on_arc(hk), base_curve);
    
    cs_.join_face(rule, true);
    //cs_.remove_vertex(v0, ->curve());
    //cs_.remove_vertex(v1);
    Halfedge_handle nh, nth;
    {
      Vertex_handle nvh= cs_.new_vertex(pt);
      nh=cs_.insert_vertex(nvh, new_extremum);
      cs_.set_curve(nh->opposite(), cs_.next_edge_on_arc(nh->opposite())->curve());
      cs_.set_curve(nh->next(), cs_.next_edge_on_arc(nh->next())->curve());
      if (nh->curve().is_inside()) nh=nh->opposite();
    }
    {
      Vertex_handle ntvh= cs_.new_vertex(Point(rule_curve, new_target->curve()));
      nth=cs_.insert_vertex(ntvh, new_extremum);
      if (nth->curve().is_inside()) nth=nth->opposite();
    }
    cs_.split_face(rule_curve, nh, nth);
  }
  // if so then remove, edge and both vertices, check edge label, add edge on other side, shoot rule
  // else roll back rule, remove rule/vertex, add rule on other side
  
}
CGAL_AOS3_END_INTERNAL_NAMESPACE
