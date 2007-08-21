#include <CGAL/Arrangement_of_spheres_3/Event_visitor.h>
#include <CGAL/Arrangement_of_spheres_3/Event_processor.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE





CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::initialize() {
  sim_= new Simulator(has_start_time_? start_time_ : tr_.min_corner()[Sweep_coordinate::index()],
		      tr_.max_corner()[Sweep_coordinate::index()]);
  
  for (CGAL_AOS3_TYPENAME Traits::Sphere_3_key_const_iterator it= tr_.sphere_3_keys_begin();
       it != tr_.sphere_3_keys_end(); ++it) {
    CGAL_AOS3_TYPENAME Traits::Event_pair ep= tr_.sphere_events(*it);
    Event_key k;
    if (ep.first >= sim_->current_time()) {
      k= sim_->new_event(ep.first,
			 CGAL_AOS3_TYPENAME Event_processor::Insert_event(j_, *it));
    }
    if (ep.second >= sim_->current_time()) {
      k= sim_->new_event(ep.second,
			 CGAL_AOS3_TYPENAME Event_processor::Remove_event(j_, *it));
    }
    if (k != Event_key() && k != sim_->null_event()) {
      free_events_.insert(k);
    }
  }
}

CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::set_traits(Traits &tr) {
  tr_=tr;
  has_traits_=true;
}

CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::set_cross_section(Combinatorial_cross_section CGAL_AOS3_TARG &cs) {
  CGAL_assertion(has_traits_);
  j_= new Event_processor CGAL_AOS3_TARG(cs, tr_);
  initialize();
}

CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::on_new_edge(Halfedge_handle h){
  if (!h->curve().key().is_input()) return;
  std::cout << "processing edge " << h->opposite()->vertex()->point()
	    << "--" << h->curve() << "--" << h->vertex()->point() << std::endl;
  
  CGAL_precondition(h->event() == Event_key());
  new_event(h);
  handle_edge_face(h, h->face());
  handle_edge_face(h, h->opposite()->face());

  process_triple(h);
  process_triple(h->opposite());
}
  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG:: on_delete_edge(Halfedge_handle h) {
  Event_key ek= h->event();
  if (ek==Event_key()) return;
  Event_key nk= sim_->null_event();
  if (ek != nk) {
    sim_->delete_event(h->event());
  }
  set_event(h, Event_key());
}
  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::on_change_edge(Halfedge_handle h) {
  if (!h->curve().key().is_input()) return;
  on_delete_edge(h);
  new_event(h);
  process_triple(h);
  process_triple(h->opposite());
}
  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::on_merge_faces(Halfedge_handle h) {
  if (!h->curve().key().is_input()) return;
  Halfedge_handle c=h;
  do {
    handle_edge_face(c, h->opposite()->face());
    c=c->next();
  } while (c != h);
} 
  

  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG:: audit(Halfedge_const_handle h) const {
  if (h->curve().key().is_input()) {
    if (!should_have_certificate(h)) {
      CGAL_assertion(h->event() == Event_key());
    } else {
      CGAL_assertion(h->event() != Event_key());
      {
	Halfedge_const_handle c=h;
	do {
	  CGAL_assertion(c->curve().key() == h->curve().key()
			 || !c->curve().key().is_input()
			 || !h->curve().key().is_input()
			 || pairs_.find(UPair(c->curve().key(),
					      h->curve().key())) != pairs_.end());
	  c=c->next();
	} while (c != h);
      }
      {
	Halfedge_const_handle c=h->opposite();
	do {
	  CGAL_assertion(c->curve().key() == h->curve().key()
			 || !c->curve().key().is_input()
			 || !h->curve().key().is_input()
			 || pairs_.find(UPair(c->curve().key(),
					      h->curve().key())) != pairs_.end());
	  c=c->next();
	} while (c != h->opposite());
      }
    
      
      if (h->next()->next()->next() ==h
	  && h->curve().is_arc()
	  && h->next()->curve().is_arc()
	  && h->next()->next()->curve().is_arc()) {
	Sphere_3_key a= h->curve().key();
	Sphere_3_key b= h->next()->curve().key();
	Sphere_3_key c= h->next()->next()->curve().key();
	if (a != b && b != c && a != c
	    && a.is_input() && b.is_input() && c.is_input()) {
	  CGAL_assertion(triples_.find(UTriple(a,b,c)) != triples_.end());
	}
      }
    }
  }
}

CGAL_AOS3_TEMPLATE
bool Event_visitor CGAL_AOS3_TARG::should_have_certificate(Halfedge_const_handle h) const {
   if (h->curve().is_arc()) {
     if ((h->vertex()->point().is_sphere_rule()
	  && h->opposite()->vertex()->point().is_sphere_rule()
	  && h->vertex()->point().rule_coordinate()
	  != h->opposite()->vertex()->point().rule_coordinate()
	  && !(h->vertex()->point().is_sphere_extremum() 
	       && h->opposite()->vertex()->point().is_sphere_extremum()))
	 || (h->vertex()->point().is_sphere_sphere() 
	     && h->opposite()->vertex()->point().is_sphere_rule()
	     && !h->opposite()->vertex()->point().is_sphere_extremum()
	     && h->opposite()->vertex()->point().rule_key() 
	     != h->vertex()->point().other_key(h->curve().key()))
	 || (h->opposite()->vertex()->point().is_sphere_sphere() 
	     && h->vertex()->point().is_sphere_rule()
	     && !h->vertex()->point().is_sphere_extremum()
	     && h->vertex()->point().rule_key() 
	     != h->opposite()->vertex()->point().other_key(h->curve().key()))){
       return true;
     } else {
       return false;
     }
   } else {
     if ( h->vertex()->point().is_rule_rule()
	 && h->opposite()->vertex()->point().is_rule_rule()
	 || (!h->next()->curve().key().is_input()
	     || !h->prev()->curve().key().is_input())) {
       return false;
     } else {
       return true;
     }
   }
}


CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::new_event(Halfedge_handle h) {
  if (!should_have_certificate(h)) {
    set_event(h, Event_key());
    std::cout << "Skipping " << h->opposite()->vertex()->point() << "--" 
	      << h->curve() << "--" << h->vertex()->point() << std::endl;
    return;
  } else {
    std::cout << "Making event for " << h->opposite()->vertex()->point() << "--" 
	      << h->curve() << "--" << h->vertex()->point() << std::endl;
  
  }
  Event_key ek=sim_->null_event();

  if (h->curve().is_arc()) {
    if (h->vertex()->point().is_sphere_rule()
	       && h->opposite()->vertex()->point().is_sphere_sphere()
	       || h->opposite()->vertex()->point().is_sphere_rule()
	       && h->vertex()->point().is_sphere_sphere()) {
      // AAR
      // circle_cross_rule_events
      Vertex_handle aav= h->vertex();
      Vertex_handle arv= h->opposite()->vertex();
      if (!aav->point().is_sphere_sphere()){
	std::swap(aav, arv);
      }
      Event_pair ep= tr_.circle_cross_rule_events(aav->point().sphere_key(0),
						  aav->point().sphere_key(1),
						  arv->point().rule_key(),
						  arv->point().rule_coordinate());
      if (ep.first.is_valid()) {
	if (ep.first >= sim_->current_time()) {
	  ek= sim_->new_event(ep.first, CGAL_AOS3_TYPENAME Event_processor::AAR_event(j_, h));
	} else if (ep.second >= sim_->current_time()) {
	  ek= sim_->new_event(ep.second, CGAL_AOS3_TYPENAME Event_processor::AAR_event(j_, h));
	} 
      } else {

      }
	
    } else if (h->vertex()->point().is_sphere_rule()
	       && h->opposite()->vertex()->point().is_sphere_rule()
	       && h->vertex()->point().rule_coordinate() 
	       != h->opposite()->vertex()->point().rule_coordinate()) {
      // RAR
      // sphere_intersect_rule_rule_events
      Sphere_3_key rks[2];
      rks[project(h->vertex()->point().rule_coordinate())]= h->vertex()->point().rule_key();
      rks[project(h->opposite()->vertex()->point().rule_coordinate())]= h->opposite()->vertex()->point().rule_key();
      Event_pair ep= tr_.sphere_intersect_rule_rule_events(h->vertex()->point().sphere_key(),
							   rks[0], rks[1]);
      if (ep.first.is_valid()) {
	if (ep.first > sim_->current_time()) {
	  ek= sim_->new_event(ep.first, CGAL_AOS3_TYPENAME Event_processor::RAR_event(j_, h));
	} else if (ep.second > sim_->current_time()) {
	  ek= sim_->new_event(ep.second, CGAL_AOS3_TYPENAME Event_processor::RAR_event(j_, h));
	} 
      }
    }
  } else {
    if (h->vertex()->point().is_sphere_rule()
	&& h->opposite()->vertex()->point().is_rule_rule()
	|| h->opposite()->vertex()->point().is_sphere_rule()
	&& h->vertex()->point().is_rule_rule()) {
      // ARR
      // sphere_intersect_rule_rule_events
      Vertex_handle arv= h->vertex();
      Vertex_handle rrv= h->opposite()->vertex();
      if (!rrv->point().is_rule_rule()){
	std::swap(rrv, arv);
      }
      Sphere_3_key rks[2];
      rks[0]= rrv->point().rule_key(plane_coordinate(0));
      rks[1]= rrv->point().rule_key(plane_coordinate(1));
   
      Event_pair ep= tr_.sphere_intersect_rule_rule_events(arv->point().sphere_key(),
							   rks[0], rks[1]);
      if (ep.first.is_valid()) {
	if (ep.first > sim_->current_time()) {
	  ek= sim_->new_event(ep.first, CGAL_AOS3_TYPENAME Event_processor::ARR_event(j_, h));
	} else if (ep.second > sim_->current_time()) {
	  ek= sim_->new_event(ep.second, CGAL_AOS3_TYPENAME Event_processor::ARR_event(j_, h));
	} else {
	  std::cout << "Both events passed: " << ep.first << " " << ep.second << std::endl;
	}
      } else {
	std::cout << "No valid event." << std::endl;
      }
    }
  }
  
  set_event(h, ek);
}
    
  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::handle_edge_face(Halfedge_handle h, Face_handle f) {
  Halfedge_handle c=f->halfedge();
  do {
    process_pair(c->curve().key(), h->curve().key());
    c=c->next();
  } while (c != f->halfedge());
}


  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::set_event(Halfedge_handle h, Event_key k) {
  h->set_event(k);
  h->opposite()->set_event(k);
}
  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::process_pair(Sphere_3_key a, 
						Sphere_3_key b) {
  if ( a != b &&  a.is_input() && b.is_input() 
       && pairs_.find(UPair(a,b))== pairs_.end()) {
    pairs_.insert(UPair(a,b));
    
    Event_pair ep= tr_.intersection_2_events(a,b);  
    if (ep.first.is_valid()) {
      CGAL_assertion(ep.first <= ep.second);
	
      if (ep.first >= sim_->current_time()) {
	Event_key k= sim_->new_event(ep.first, CGAL_AOS3_TYPENAME Event_processor::Intersect_event(j_, a,b));
	if (k != Event_key() && k != sim_->null_event()) {
	  free_events_.insert(k);
	}
      }
      if (ep.second >= sim_->current_time()) {
	Event_key k =sim_->new_event(ep.second, CGAL_AOS3_TYPENAME Event_processor::Unintersect_event(j_, a,b));
	if (k != Event_key() && k != sim_->null_event()) {
	  free_events_.insert(k);
	}
      }
	
	
      // now handle the extremum intersections
      for (unsigned int k=0; k< 2; ++k) {
	for (unsigned int i=0; i< 2; ++i){
	  Event_pair ep= tr_.sphere_intersect_extremum_events(a,
							      plane_coordinate(i),
							      b);  
	  if (ep.first.is_valid()) {
	    CGAL_assertion(ep.first <= ep.second);
	    if (ep.first >= sim_->current_time()) {
	      Event_key k= sim_->new_event(ep.first, CGAL_AOS3_TYPENAME Event_processor::AAE_event(j_, a,b, i));
	      if (k != Event_key() && k != sim_->null_event()) {
		free_events_.insert(k);
	      }
	    }
	    if (ep.second >= sim_->current_time()) {
	      Event_key k =sim_->new_event(ep.second, CGAL_AOS3_TYPENAME Event_processor::AAE_event(j_, a,b, i));
	      if (k != Event_key() && k != sim_->null_event()) {
		free_events_.insert(k);
	      }
	    }
	  }
	}
	std::swap(a,b);
      }
    }
  }
}


CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::process_triple(Halfedge_handle h) {
  if (!h->curve().is_arc() || !h->next()->curve().is_arc() || !h->next()->next()->curve().is_arc()) {
    if (h->curve().is_arc()
	&& h->next()->curve().is_arc()
	&& h->next()->next()->curve().is_arc()) {
      Sphere_3_key a= h->curve().key();
      Sphere_3_key b= h->next()->curve().key();
      Sphere_3_key c= h->next()->next()->curve().key();
     
      if (a !=b && b!= c && a != c 
	  && a.is_input() && b.is_input() && c.is_input() 
	  && triples_.find(UTriple(a,b,c)) == triples_.end()) {
	triples_.insert(UTriple(a,b,c));
	Event_pair ep= tr_.intersection_3_events(a,b,c);  
	if (ep.first.is_valid()) {
	  CGAL_assertion(ep.first <= ep.second);
	  if (ep.first >= sim_->current_time()) {
	    Event_key k= sim_->new_event(ep.first, CGAL_AOS3_TYPENAME Event_processor::Intersect_3_event(j_, a,b,c));
	    if (k != Event_key() && k != sim_->null_event()) {
	      free_events_.insert(k);
	    }
	  }
	  if (ep.second >= sim_->current_time()) {
	    Event_key k= sim_->new_event(ep.second, CGAL_AOS3_TYPENAME Event_processor::Unintersect_3_event(j_, a,b,c));
	    if (k != Event_key() && k != sim_->null_event()) {
	      free_events_.insert(k);
	    }
	  }
	}
      }
    }
  }
}

CGAL_AOS3_END_INTERNAL_NAMESPACE
