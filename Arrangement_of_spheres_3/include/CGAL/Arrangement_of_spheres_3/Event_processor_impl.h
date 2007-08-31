#include <CGAL/Arrangement_of_spheres_3/Event_processor.h>
#include <algorithm>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

#define CGAL_CATCH_DEGENERACY int
//Degeneracy

CGAL_AOS3_TEMPLATE
Event_processor CGAL_AOS3_TARG::Event_processor(CCS &cs, Traits &tr): cs_(cs), tr_(tr){
  
}

CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::handle_degeneracy() {  
  CGAL_AOS3_TYPENAME Traits::Event_point_3 ep= cs_.visitor().simulator()->current_time();
  std::cout << "Degeneracy at " << ep << std::endl;
  std::vector<Vertex_handle> vertices;
  ICSL ics(tr_, cs_);
  for (CGAL_AOS3_TYPENAME CCS::Vertex_const_iterator it= cs_.vertices_begin();
       it != cs_.vertices_end(); ++it) {
    if (ics.equal_points(it, ep)) {
      std::cout << it->point() << std::endl;
    }
  }
  //CGAL_assertion(0);
}


CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::check_degeneracy() {  
  if (cs_.visitor().simulator()->current_time() 
      == cs_.visitor().simulator()->next_event_time()) {
    std::cout << "Degeneracy at " << cs_.visitor().simulator()->current_time() << std::endl;
    throw Degeneracy();
  } else {
    std::cout << "Current time is " << cs_.visitor().simulator()->current_time()  
	      << " and next event time is " << cs_.visitor().simulator()->next_event_time() << std::endl;
  }
}


CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::insert(Sphere_3_key k) {  
  cs_.audit();
  try {
    check_degeneracy();
    ICSL icsl(tr_, cs_);
    ICSI icsi(tr_, cs_);
    try {
      CGAL_AOS3_TYPENAME CCS::Face_handle f= icsl.locate(k);
      //slice.new_marked_face(f);
      icsi.insert(k, f);
    } catch (CGAL_AOS3_TYPENAME ICSL::On_edge_exception e) {
      icsi.insert(k, e.halfedge_handle());
    } catch (CGAL_AOS3_TYPENAME ICSL::On_vertex_exception v) {
      icsi.insert(k, v.vertex_handle());
    }
  } catch (CGAL_CATCH_DEGENERACY) {
    std::cout << "Degeneracy" << std::endl;
    handle_degeneracy();
  }
  cs_.audit();
}

CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::remove(Sphere_3_key k) {
  cs_.audit();
  try {
    // just try to process it anyway and see if I can
    check_degeneracy();
    ICSD icsr(tr_, cs_);
    icsr.remove_sphere(k);
  } catch (CGAL_CATCH_DEGENERACY) {
    std::cout << "Degeneracy" << std::endl;
    handle_degeneracy();
  }
  cs_.audit();
}

CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::intersect(Sphere_3_key k, Sphere_3_key l) {
  cs_.audit();
  try {
    check_degeneracy();

    if (tr_.compare_sphere_centers_c(k,l, plane_coordinate(0)) == CGAL::EQUAL
	|| tr_.compare_sphere_centers_c(k,l, plane_coordinate(1)) == CGAL::EQUAL) {
      std::cout << "Centers line up" << std::endl;
      throw Degeneracy();
    }
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
      ICSL ics(tr_, cs_);
    
      f=ics.locate(border_both.begin(), border_both.end(),
		   cs_.visitor().simulator()->current_time());
      //f= ics.locate(border_both.begin(), border_both.end());
    } else {
      f= border_both.back();
    }

    // have problem with more than one vertex on face

    std::cout << "Intersecting on face ";
    cs_.write(f, std::cout) << std::endl;

    Halfedge_handle hk, hl;
    
    {
      
      std::vector<Halfedge_handle> hk_cand, hl_cand;
      
      Halfedge_handle hc=f->halfedge();
      do {
	if (hc->curve().is_arc()) {
	  if (hc->curve().key() == k) {
	    hk_cand.push_back(hc);
	  } else if (hc->curve().key() == l) {
	    hl_cand.push_back(hc);
	  }
	}
	hc=hc->next();
      } while (hc != f->halfedge());

      if (hk_cand.size() > 1 || hl_cand.size() >1) {
	ICSR icsr(tr_, cs_);
	for (int i=0; i< 2; ++i) {
	  try {
	    Halfedge_handle h= icsr.shoot_rule(cs_.visitor().simulator()->current_time(),
					       f,
					       Rule_direction(i*2));
	    CGAL_assertion(h->curve().is_arc());
	    if (h->curve.key()== k) {
	      hk=h;
	    } else {
	      CGAL_assertion(h->curve().key() ==l);
	      hl=h;
	    }
	  } catch (CGAL_AOS3_TYPENAME ICSR::On_vertex_exception) {
	    throw Degeneracy_exception();
	  }
	}
      }
    }
    std::cout << "Intersection edges are ";
    cs_.write(hk, std::cout) << " and ";
    cs_.write(hl, std::cout) << std::endl;

    

    cs_.intersect_arcs(hl, hk);
  } catch (CGAL_CATCH_DEGENERACY) {
    std::cout << "Degeneracy" << std::endl;
    handle_degeneracy();
  }
  cs_.audit();
}







CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::unintersect(Sphere_3_key k, Sphere_3_key l) {
  cs_.audit();
  try {
    check_degeneracy();
    if (tr_.compare_sphere_centers_c(k,l, plane_coordinate(0)) == CGAL::EQUAL
	|| tr_.compare_sphere_centers_c(k,l, plane_coordinate(1)) == CGAL::EQUAL) {
      std::cout << "Centers line up" << std::endl;
      throw Degeneracy();
    }
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
      hk= cs_.next_edge_on_circle(hk);
    } while (hk != khs);
    if (hl == Halfedge_handle() || hk== Halfedge_handle()) {
      throw Degeneracy();
    }

    {
      ICSL icsl(tr_, cs_);
      CGAL_assertion(icsl.equal_points(hk->vertex(), cs_.visitor().simulator()->current_time()));
      CGAL_assertion(icsl.equal_points(hl->vertex(), cs_.visitor().simulator()->current_time()));
      CGAL_assertion(hk->vertex()== hl->opposite()->vertex());
      CGAL_assertion(hl->vertex()== hk->opposite()->vertex());
    }

    // find shared face, f, and edges hk and hl
    cs_.unintersect_arcs(hl, hk);
  } catch (Degeneracy) {
    std::cout << "Degeneracy" << std::endl;
    handle_degeneracy();
  }
  cs_.audit();
}








CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::intersect(Sphere_3_key k, Sphere_3_key l, Sphere_3_key m) {
  cs_.audit();
  ICSL icsl(tr_, cs_);

  try {
    check_degeneracy();
    
    std::vector<Halfedge_handle> cands;
    {
      Halfedge_handle hk=cs_.a_halfedge(k), khs=hk;
      
      do {
	if (hk->next()->next()->next() == hk) {
	  if ((hk->next()->curve().key() == l 
	       || hk->next()->next()->curve().key() == l)
	    && (hk->next()->curve().key() == m
		|| hk->next()->next()->curve().key() == m)){
	    cands.push_back(hk);
	  }
	} else if (hk->opposite()->next()->next()->next()->opposite()==hk){
	  if ((hk->opposite()->next()->curve().key() == l 
	       || hk->opposite()->next()->next()->curve().key() == l)
	      && (hk->opposite()->next()->curve().key() == m
		  || hk->opposite()->next()->next()->curve().key() == m)){
	    cands.push_back(hk->opposite());
	  }
	} 
	hk= cs_.next_edge_on_circle(hk);
      } while (hk != khs);
      
      if (cands.empty()) {
	throw Degeneracy();
      }
    }

    Halfedge_handle hk;
    if (cands.size() >1) {
      for (unsigned int i=0; i< cands.size(); ++i){
	if (icsl.equal_points(cands[i]->vertex(), cs_.visitor().simulator()->current_time())
	    && icsl.equal_points(cands[i]->next()->vertex(), cs_.visitor().simulator()->current_time())
	    &&  icsl.equal_points(cands[i]->next()->next()->vertex(), cs_.visitor().simulator()->current_time())){
	  CGAL_assertion(hk== Halfedge_handle());
	  hk= cands[i];
#ifdef NDEBUG
	  break;
#endif
	}
      }
    } else {
      hk=cands.back();
      cands.clear();
    }

    {
      //ICSL icsl(tr_, cs_);
      CGAL_assertion(icsl.equal_points(hk->vertex(), cs_.visitor().simulator()->current_time()));
      CGAL_assertion(icsl.equal_points(hk->next()->vertex(), cs_.visitor().simulator()->current_time()));
      CGAL_assertion(icsl.equal_points(hk->next()->next()->vertex(), cs_.visitor().simulator()->current_time()));
      //CGAL_assertion(hk->vertex()== hl->opposite()->vertex());
      //CGAL_assertion(hl->vertex()== hk->opposite()->vertex());
      CGAL_assertion(hk->curve().key() ==k);
      CGAL_assertion(hk->next()->curve().key() == l 
		     || hk->next()->next()->curve().key() == l);
      CGAL_assertion(hk->next()->curve().key() == m
		     || hk->next()->next()->curve().key() == m);
    }

    std::cout << "Triple edges are \n";
    cs_.write(hk, std::cout) << std::endl;
    cs_.write(hk->next(), std::cout) << std::endl;
    cs_.write(hk->next()->next(), std::cout) << std::endl;

    Halfedge_handle hkn= cs_.next_edge_on_circle(hk)->opposite();
    Halfedge_handle hkp= cs_.next_edge_on_circle(hk->opposite())->opposite();
  
    cs_.write(hkn, std::cout) << " is hkn\n";
    cs_.write(hkp, std::cout) << " is hkp\n";
    
    Halfedge_handle hknt= hkn->opposite()->prev()->prev();
    Halfedge_handle hkpt= hkp->next()->next();

    cs_.write(hknt, std::cout) << " is hknt\n";
    cs_.write(hkpt, std::cout) << " is hkpt\n";
  
    Point pn= hkp->vertex()->point(), pp= hkn->vertex()->point();
  
    std::cout << "Points are " << pn << " and " << pp << std::endl;
    Vertex_handle nvhn= cs_.new_vertex(pn);
    Vertex_handle nvhp= cs_.new_vertex(pp);

    cs_.insert_vertex(nvhn, hknt);
    cs_.insert_vertex(nvhp, hkpt);
    cs_.move_target(hkp, nvhp, true);
    cs_.move_target(hkn, nvhn, true);


    cs_.split_face(hk->curve(),  
		   hkp->next()->opposite(),
		   hkn->next()->opposite());
    cs_.join_face(hk, true);
  } catch (CGAL_CATCH_DEGENERACY) {
    std::cout << "Degeneracy" << std::endl;
    handle_degeneracy();
  }
  cs_.audit();
}







CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::unintersect(Sphere_3_key k, Sphere_3_key l, Sphere_3_key m) {
  intersect(k,l,m);
}







CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::process_aar(Halfedge_handle h) {
  cs_.audit();
  try {
    check_degeneracy();
    clear_event(h);

    CGAL_assertion(h->vertex()->point().is_sphere_rule());
    CGAL_assertion(h->opposite()->vertex()->point().is_sphere_sphere());
    tr_.advance_circle_cross_rule_event(h->opposite()->vertex()->point().sphere_key(0),
					h->opposite()->vertex()->point().sphere_key(1),
					h->vertex()->point().rule_key(),
					h->vertex()->point().rule_constant_coordinate());
    std::cout << "advance AAR of " << h->opposite()->vertex()->point().sphere_key(0) << " " 
	      << h->opposite()->vertex()->point().sphere_key(1) << " "
	      << h->vertex()->point().rule_key() << " "
	      << h->vertex()->point().rule_constant_coordinate() << std::endl;

    {
      ICSL icsl(tr_, cs_);
      CGAL_assertion(icsl.equal_points(h->vertex(), cs_.visitor().simulator()->current_time()));
      CGAL_assertion(icsl.equal_points(h->opposite()->vertex(), cs_.visitor().simulator()->current_time()));
    }

    Halfedge_handle rule= cs_.cross_edge(h)->opposite();
    CGAL_assertion(rule->vertex() ==h->vertex());
    CGAL_assertion(rule->curve().is_rule());
    Halfedge_handle target;
    if (h->next() == rule->opposite()) {
      target= h->prev();
    } else {
      target=h->opposite()->next();
    }
    
    Vertex_handle nv= cs_.new_vertex(Combinatorial_vertex(target->curve(),
							  rule->curve()));
    std::cout << "New vertex of " << nv->point() << std::endl;
    cs_.insert_vertex(nv, target);
    std::cout << "Inserted" << std::endl;
    cs_.move_target(rule, nv, true);
    cs_.audit();
  } catch (CGAL_CATCH_DEGENERACY) {
    std::cout << "Degeneracy" << std::endl;
    handle_degeneracy();
  }
}







CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::process_rar(Halfedge_handle h) {
  cs_.audit();
  try {
    check_degeneracy();
    clear_event(h);

    {
      ICSL icsl(tr_, cs_);
      CGAL_assertion(icsl.equal_points(h->vertex(), cs_.visitor().simulator()->current_time()));
      CGAL_assertion(icsl.equal_points(h->opposite()->vertex(), cs_.visitor().simulator()->current_time()));
    }

    Sphere_3_key rks[2];
    rks[project(h->vertex()->point().rule_constant_coordinate())]= h->vertex()->point().rule_key();
    rks[project(h->opposite()->vertex()->point().rule_constant_coordinate())]
      = h->opposite()->vertex()->point().rule_key();
    
    tr_.advance_sphere_intersect_rule_rule_event(h->vertex()->point().sphere_key(), rks[0], rks[1]);
						 
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
  } catch (CGAL_CATCH_DEGENERACY) {
    std::cout << "Degeneracy" << std::endl;
    handle_degeneracy();
  }
  cs_.audit();
}


CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::clear_event(Halfedge_handle h) {
  h->set_event(cs_.visitor().simulator()->null_event());
  h->opposite()->set_event(cs_.visitor().simulator()->null_event());
}
  

CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::process_arr(Halfedge_handle h) {
  cs_.audit();
  try {
    check_degeneracy();
    clear_event(h);
    
    {
      CGAL_precondition(h->vertex()->point().is_rule_rule());
      ICSL icsl(tr_, cs_);
      CGAL_assertion(icsl.equal_points(h->vertex(), cs_.visitor().simulator()->current_time()));
      CGAL_assertion(icsl.equal_points(h->opposite()->vertex(), cs_.visitor().simulator()->current_time()));
    }
    
    
    Sphere_3_key rks[2];
    rks[0]= h->vertex()->point().rule_key(plane_coordinate(0));
    rks[1]= h->vertex()->point().rule_key(plane_coordinate(1));
   
    tr_.advance_sphere_intersect_rule_rule_event(h->opposite()->vertex()->point().sphere_key(),
						 rks[0], rks[1]);


    if (!h->opposite()->vertex()->point().is_rule_rule()) h= h->opposite();
    CGAL_assertion(h->opposite()->vertex()->point().is_rule_rule());
    Halfedge_handle xr= cs_.cross_edge(h->opposite())->opposite(); // inward pointing
    //CGAL_assertion(xr->vertex() == h->vertex());
    std::cout << "XR is ";
    cs_.write(xr, std::cout) << std::endl;

    ICSR ics(tr_, cs_);
    ics.roll_back_rule(cs_.visitor().simulator()->current_time(), xr);
    ics.roll_back_rule(cs_.visitor().simulator()->current_time(), xr->opposite());
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
    }
    if (cs_.is_redundant(xr)) {
      std::cout << "Joining redundant faces on edge ";
      cs_.write(xr, std::cout) << std::endl;
      cs_.join_face(xr, true);
    } else {
      do {
	Halfedge_handle target_h;
	if (h->opposite()->face() == xr->face() || h->opposite()->face() == xr->opposite()->face()) {
	  target_h= h->opposite()->prev();
	} else {
	  target_h=h->next();
	}
	std::cout << "target is ";
	cs_.write(target_h, std::cout) << std::endl;
	
	Point pt(xr->curve(), target_h->curve());
	Vertex_handle nv= cs_.new_vertex(pt);
	cs_.insert_vertex(nv, target_h);
	cs_.move_target(xr, nv, true);
      } while ((xr=nxr) != Halfedge_handle()); // do other side if needed
    }
  } catch (CGAL_CATCH_DEGENERACY) {
    std::cout << "Degeneracy" << std::endl;
    handle_degeneracy();
  }
  cs_.audit();
  // could fail if there is a degeneracy

  // now the cross rule should be a t pointing on direction;
  
 
}

CGAL_AOS3_TEMPLATE
void Event_processor CGAL_AOS3_TARG::process_aae(Sphere_3_key k, Sphere_3_key l, int i) {
  cs_.audit();
  try {
    check_degeneracy();
    int index=i+1;
    bool larger=tr_.compare_point_to_rule_c(cs_.visitor().simulator()->current_time(), k,
				      plane_coordinate(1-i))== CGAL::LARGER;
    if (!larger && i==0 || larger && i==1) {
      index=(index+2)%4;
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
			       << " " <<  (hk->vertex()->point().is_sphere_extremum()? 
					   hk->vertex()->point().sphere_extremum_index().index()
					   : -1) << std::endl;
      hk= cs_.next_edge_on_circle(hk);
    }
    Halfedge_handle rule= hk->opposite()->prev()->opposite(); // outward rule
 

    std::cout << index << " hk is ";
    cs_.write(hk, std::cout) << " and rule is ";
    cs_.write(rule, std::cout) << std::endl;

    if (tr_.compare_sphere_centers_c(k,l, plane_coordinate(i)) == CGAL::EQUAL) {
      std::cout << "Centers line up" << std::endl;
      throw Degeneracy();
    }


    {
      ICSR icsr(tr_, cs_);
      ICSL icsl(tr_, cs_);
      icsr.roll_back_rule(cs_.visitor().simulator()->current_time(), rule);
      std::cout << "hk is now ";
      cs_.write(hk, std::cout) << " and rule is ";
      cs_.write(rule, std::cout) << std::endl;

      Halfedge_handle middle_edge; // in middle point to rule
      Face_handle new_face;
      Curve rule_curve=rule->curve();
      Point extremum= hk->vertex()->point();
    
      CGAL_precondition(icsl.equal_points(hk->vertex(), cs_.visitor().simulator()->current_time()));

      if (icsl.equal_points(hk->opposite()->vertex(), cs_.visitor().simulator()->current_time())) {
	if (icsl.equal_points(cs_.next_edge_on_circle(hk)->vertex(), cs_.visitor().simulator()->current_time())) {
	  Degeneracy d;
	  d.new_vertex(hk->opposite()->vertex());
	  d.new_vertex(cs_.next_edge_on_circle(hk)->vertex());
	  throw d; 
	}
	middle_edge= hk;
      } else {
	middle_edge=cs_.next_edge_on_circle(hk)->opposite();
	std::cout << "Middle edge is ";
	cs_.write(middle_edge, std::cout) << std::endl;
	CGAL_assertion(icsl.equal_points(middle_edge->opposite()->vertex(), 
					 cs_.visitor().simulator()->current_time()));
	//middle_edge=cs_.next_edge_on_circle(hk)->opposite();
      }
      std::cout << "Middle is ";
      cs_.write(middle_edge, std::cout) << std::endl;

      
      Halfedge_handle new_extremum_location= cs_.next_edge_on_circle(middle_edge->opposite());
      Face_handle search_face;
      if (new_extremum_location->curve().is_inside()) search_face= new_extremum_location->opposite()->face();
      else search_face=  new_extremum_location->face();

      std::cout << "New extremum is ";
      cs_.write(new_extremum_location, std::cout) << std::endl;

      cs_.set_curve(middle_edge, cs_.next_edge_on_circle(middle_edge)->curve());
      cs_.join_face(rule, true);
   
      Vertex_handle nvh= cs_.new_vertex(extremum);
      cs_.insert_vertex(nvh, new_extremum_location);
      Halfedge_handle source_h= cs_.find_halfedge(nvh, search_face);


      Halfedge_handle hd= icsr.find_rule_vertex(cs_.visitor().simulator()->current_time(),
					       search_face,
					       rule_curve);
      cs_.split_face(rule_curve, source_h, hd);
    }

    // if so then remove, edge and both vertices, check edge label, add edge on other side, shoot rule
    // else roll back rule, remove rule/vertex, add rule on other side
     cs_.audit();
  } catch (CGAL_CATCH_DEGENERACY) {
    std::cout << "Degeneracy" << std::endl;
    handle_degeneracy();
  }
}
CGAL_AOS3_END_INTERNAL_NAMESPACE
