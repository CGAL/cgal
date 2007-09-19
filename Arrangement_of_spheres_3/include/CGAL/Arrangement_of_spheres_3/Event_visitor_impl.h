#include <CGAL/Arrangement_of_spheres_3/Event_visitor.h>
#include <CGAL/Arrangement_of_spheres_3/Event_processor.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE





CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::initialize() {
  sim_= new Simulator(has_start_time_? start_time_ : tr_.min_corner()[Sweep_coordinate::index()],
		      tr_.max_corner()[Sweep_coordinate::index()]+10);
  
  for (CGAL_AOS3_TYPENAME Traits::Sphere_3_key_const_iterator it= tr_.sphere_3_keys_begin();
       it != tr_.sphere_3_keys_end(); ++it) {
    CGAL_AOS3_TYPENAME Traits::Event_pair ep= tr_.sphere_events(*it);
    Event_key k;
    if (ep.first >= sim_->current_time()) {
      k= sim_->new_event(ep.first,
			 CGAL_AOS3_TYPENAME Event_processor CGAL_AOS3_TARG::I_event(j_, *it));
    }
    if (ep.second >= sim_->current_time()) {
      k= sim_->new_event(ep.second,
			 CGAL_AOS3_TYPENAME Event_processor CGAL_AOS3_TARG::R_event(j_, *it));
    }
    if (k != Event_key() && k != sim_->null_event()) {
      free_events_.insert(k);
    }
  }

  {
 
    std::vector<Box> boxes;
    for (CGAL_AOS3_TYPENAME Traits::Sphere_3_key_const_iterator it
	   = tr_.sphere_3_keys_begin(); it != tr_.sphere_3_keys_end(); ++it) {
      boxes.push_back(Box(tr_.sphere_3(*it).bbox(), std::make_pair(*it, *it)));
    }

    CGAL::box_self_intersection_d(boxes.begin(), boxes.end(), 
				  Report_pairs(this));

    //std::copy(intersections_.begin(), intersections_.end(), std::back_inserter(boxes));
    CGAL::box_intersection_d(boxes.begin(), boxes.end(),
			     intersections_.begin(), intersections_.end(),
			     Report_triples(this));
    intersections_.clear();
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
  CGAL_LOG(Log::LOTS, "processing edge " << h->opposite()->vertex()->point()
	   << "--" << h->curve() << "--" << h->vertex()->point() << std::endl);
  
  CGAL_precondition(h->event() == Event_key());
  new_event(h);
  handle_edge_face(h, h->face());
  handle_edge_face(h, h->opposite()->face());

  //process_triple(h);
  //process_triple(h->opposite());
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
  //process_triple(h);
  //process_triple(h->opposite());
}
  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::on_merge_faces(Halfedge_handle ) {
  //if (!h->curve().key().is_input()) return;
  /*Halfedge_handle c=h;
  do {
    if (c->curve().key().is_input()) {
      handle_edge_face(c, h->opposite()->face());
    }
    c=c->next();
    } while (c != h);*/
} 
  

  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG:: audit(Halfedge_const_handle h) const {
  if (h->curve().key().is_input()) {
    if (!should_have_certificate(h)) {
      CGAL_assertion(h->event() == Event_key());
    } else {
      /*CGAL_assertion(h->event() != Event_key());
      {
	Halfedge_const_handle c=h;
	do {
	  CGAL_assertion(c->curve().key() == h->curve().key()
			 || !c->curve().key().is_input()
			 || !h->curve().key().is_input()
			 || (!tr_.intersects(c->curve().key(),h->curve().key())
			     || pairs_.find(Sphere_key_upair(c->curve().key(),
							     h->curve().key())) != pairs_.end()));
	  c=c->next();
	} while (c != h);
      }
      {
	Halfedge_const_handle c=h->opposite();
	do {
	  CGAL_assertion(c->curve().key() == h->curve().key()
			 || !c->curve().key().is_input()
			 || !h->curve().key().is_input()
			 || (!tr_.intersects(c->curve().key(),h->curve().key())
			     || pairs_.find(Sphere_key_upair(c->curve().key(),
							     h->curve().key())) != pairs_.end()));
	  c=c->next();
	} while (c != h->opposite());
	}*/
    
      
      if (h->next()->next()->next() ==h
	  && h->curve().is_arc()
	  && h->next()->curve().is_arc()
	  && h->next()->next()->curve().is_arc()) {
	Sphere_3_key a= h->curve().key();
	Sphere_3_key b= h->next()->curve().key();
	Sphere_3_key c= h->next()->next()->curve().key();
	/*if (a != b && b != c && a != c
	    && a.is_input() && b.is_input() && c.is_input()) {
	  CGAL_assertion(triples_.find(Sphere_key_utriple(a,b,c)) != triples_.end());
	  }*/
      }
    }
  }
}

CGAL_AOS3_TEMPLATE
bool Event_visitor CGAL_AOS3_TARG::should_have_certificate(Halfedge_const_handle h) const {
   if (h->curve().is_arc()) {
     if ((h->vertex()->point().is_sphere_rule()
	  && h->opposite()->vertex()->point().is_sphere_rule()
	  && h->vertex()->point().rule_constant_coordinate()
	  != h->opposite()->vertex()->point().rule_constant_coordinate()
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
    CGAL_LOG(Log::LOTS, "Skipping " << h->opposite()->vertex()->point() << "--" 
	     << h->curve() << "--" << h->vertex()->point() << std::endl);
    return;
  } else {
    CGAL_LOG(Log::LOTS, "Making event for " << h->opposite()->vertex()->point() << "--" 
	     << h->curve() << "--" << h->vertex()->point() << std::endl);
  
  }
  Event_key ek=sim_->null_event();

  if (h->curve().is_arc()) {
    if (h->vertex()->point().is_sphere_rule()
	       && h->opposite()->vertex()->point().is_sphere_sphere()
	       || h->opposite()->vertex()->point().is_sphere_rule()
	       && h->vertex()->point().is_sphere_sphere()) {
      // AAR
      // circle_cross_rule_events
      if (!h->vertex()->point().is_sphere_rule()) h=h->opposite();
      Vertex_handle aav= h->opposite()->vertex();
      Vertex_handle arv= h->vertex();
      /*if (!aav->point().is_sphere_sphere()){
	std::swap(aav, arv);
	}*/

      do {
	CGAL_LOG(Log::LOTS, "create AAR of " << aav->point().sphere_key(0) << " " 
		  << aav->point().sphere_key(1) << " "
		  << arv->point().rule_key() << " "
		 << arv->point().rule_constant_coordinate() << std::endl);
	CGAL_AOS3_TYPENAME Traits::Event_point_3 ep= tr_.circle_cross_rule_event(aav->point().sphere_key(0),
										 aav->point().sphere_key(1),
										 arv->point().rule_key(),
										 arv->point().rule_constant_coordinate());
	if (ep.is_valid()) {
	  //if (tr_.equal(ep, sim_->current_time())) throw CGAL_AOS3_TYPENAME Traits::Degeneracy_exception();
	  if (ep >= sim_->current_time()) {
	    CGAL_assertion(tr_.oriented_side_of_separating_plane(ep, aav->point().sphere_key(0),
								       aav->point().sphere_key(1)) != CGAL::NEGATIVE);
	    ek= sim_->new_event(ep, CGAL_AOS3_TYPENAME EP::AAR_event(j_, h));
	    break;
	  } else {
	    CGAL_LOG(Log::LOTS, "advance AAR of " << aav->point().sphere_key(0) << " " 
		      << aav->point().sphere_key(1) << " "
		      << arv->point().rule_key() << " "
		     << arv->point().rule_constant_coordinate() << std::endl);

	    tr_.advance_circle_cross_rule_event(aav->point().sphere_key(0),
						aav->point().sphere_key(1),
						arv->point().rule_key(),
						arv->point().rule_constant_coordinate());
	  }
	} else break;
      } while (true);
	
    } else if (h->vertex()->point().is_sphere_rule()
	       && h->opposite()->vertex()->point().is_sphere_rule()
	       && h->vertex()->point().rule_constant_coordinate() 
	       != h->opposite()->vertex()->point().rule_constant_coordinate()) {
      // RAR
      // sphere_intersect_rule_rule_events
      Sphere_3_key rks[2];
      rks[project(h->vertex()->point().rule_constant_coordinate())]= h->vertex()->point().rule_key();
      rks[project(h->opposite()->vertex()->point().rule_constant_coordinate())]= h->opposite()->vertex()->point().rule_key();
      do {
	CGAL_AOS3_TYPENAME Traits::Event_point_3 ep= tr_.sphere_intersect_rule_rule_event(h->vertex()->point().sphere_key(),
											  rks[0], rks[1]);
	if (ep.is_valid()) {
	  //if (tr_.equal(ep, sim_->current_time())) throw CGAL_AOS3_TYPENAME Traits::Degeneracy_exception();
	  if (ep >= sim_->current_time()) {
	    ek= sim_->new_event(ep, CGAL_AOS3_TYPENAME EP::RAR_event(j_, h));
	    break;
	  } else {
	    tr_.advance_sphere_intersect_rule_rule_event(h->vertex()->point().sphere_key(),
							 rks[0], rks[1]);
	  }
	} else break;
      } while (true);
    }
  } else {
    if (h->vertex()->point().is_sphere_rule()
	&& h->opposite()->vertex()->point().is_rule_rule()
	|| h->opposite()->vertex()->point().is_sphere_rule()
	&& h->vertex()->point().is_rule_rule()) {
      // ARR
      // sphere_intersect_rule_rule_events
      if (h->vertex()->point().is_sphere_rule()) h=h->opposite();
      Vertex_handle arv= h->opposite()->vertex();
      Vertex_handle rrv= h->vertex();
      /*if (!rrv->point().is_rule_rule()){
	std::swap(rrv, arv);
	h=h->opposite();
	}*/
      Sphere_3_key rks[2];
      rks[0]= rrv->point().rule_key(plane_coordinate(0));
      rks[1]= rrv->point().rule_key(plane_coordinate(1));
   
      do {
	CGAL_AOS3_TYPENAME Traits::Event_point_3 ep= tr_.sphere_intersect_rule_rule_event(arv->point().sphere_key(),
											  rks[0], rks[1]);
	if (ep.is_valid()) {
	  //if (tr_.equal(ep, sim_->current_time())) throw CGAL_AOS3_TYPENAME Traits::Degeneracy_exception();
	  if (ep >= sim_->current_time()) {
	    ek= sim_->new_event(ep, CGAL_AOS3_TYPENAME EP::ARR_event(j_, h));
	    break;
	  } else {
	    tr_.advance_sphere_intersect_rule_rule_event(arv->point().sphere_key(),
							  rks[0], rks[1]);
	  }
	} else {
	  break;
	}
      } while (true);
    }
  }
  
  set_event(h, ek);
}
    
  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::handle_edge_face(Halfedge_handle , Face_handle ) {
  /*Halfedge_handle c=f->halfedge();
  do {
    if (c->curve().key().is_input()) {
      process_pair(c->curve().key(), h->curve().key());
    }
    c=c->next();
    } while (c != f->halfedge());*/
}


  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::set_event(Halfedge_handle h, Event_key k) {
  h->set_event(k);
  h->opposite()->set_event(k);
}
  
CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::process_pair(const Box & a, 
						const Box & b) {
  if (tr_.intersects(a.handle().first, b.handle().first)) {
    Bbox_3 bb(std::max(a.min_coord(0), b.min_coord(0)),
	      std::max(a.min_coord(1), b.min_coord(1)),
	      std::max(a.min_coord(2), b.min_coord(2)),
	      std::min(a.max_coord(0), b.max_coord(0)),
	      std::min(a.max_coord(1), b.max_coord(1)),
	      std::min(a.max_coord(2), b.max_coord(2)));
    intersections_.push_back(Box(bb, std::make_pair(a.handle().first, b.handle().first)));
    process_pair(a.handle().first, b.handle().first);
  }
}


CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::process_pair(Sphere_3_key a, 
						Sphere_3_key b) {
  if ( a != b &&  a.is_input() && b.is_input() 
       /*&& pairs_.find(Sphere_key_upair(a,b))== pairs_.end()*/) {
    CGAL_LOG(Log::LOTS, "Processing pair " << a << " " << b << std::endl);
    // pairs_.insert(Sphere_key_upair(a,b));
    Event_pair ep= tr_.intersection_2_events(a,b);
   
    if (ep.first.is_valid() && !ep.second.is_valid()) {
      CGAL_LOG(Log::SOME, "Spheres " << a << " and " 
	       << b << " intersect on vertical circle"<< std::endl);
      Event_key k= sim_->new_event(ep.first, CGAL_AOS3_TYPENAME EP::S2_event(j_, a,b));
      if (k != Event_key() && k != sim_->null_event()) {
	free_events_.insert(k);
      }
    } else if (ep.first.is_valid()) {
	CGAL_assertion(ep.first <= ep.second);
	
	if (ep.first >= sim_->current_time()) {
	  Event_key k= sim_->new_event(ep.first, CGAL_AOS3_TYPENAME EP::I2_event(j_, a,b));
	  if (k != Event_key() && k != sim_->null_event()) {
	    free_events_.insert(k);
	  }
	}
	if (ep.second >= sim_->current_time()) {
	  Event_key k =sim_->new_event(ep.second, CGAL_AOS3_TYPENAME EP::U2_event(j_, a,b));
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
		Event_key k= sim_->new_event(ep.first, CGAL_AOS3_TYPENAME EP::AAE_event(j_, a,b, i));
		if (k != Event_key() && k != sim_->null_event()) {
		  free_events_.insert(k);
		}
	      }
	      if (ep.second >= sim_->current_time()) {
		Event_key k =sim_->new_event(ep.second, CGAL_AOS3_TYPENAME EP::AAE_event(j_, a,b, i));
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
void Event_visitor CGAL_AOS3_TARG::process_triple(Sphere_3_key a, 
						  Sphere_3_key b,
						  Sphere_3_key c) {
  if (a !=b && b!= c && a != c 
      && a.is_input() && b.is_input() && c.is_input()) {
    //triples_.insert(Sphere_key_utriple(a,b,c));
    Event_pair ep= tr_.intersection_3_events(a,b,c);  
    if (ep.first.is_valid()) {
      CGAL_assertion(ep.first <= ep.second);
      if (a>b) std::swap(a,b);
      if (b>c) std::swap(b,c);
      if (a>b) std::swap(a,b);
      
      if (ep.first >= sim_->current_time()) {
	Event_key k= sim_->new_event(ep.first, CGAL_AOS3_TYPENAME EP::I3_event(j_, a,b,c));
	if (k != Event_key() && k != sim_->null_event()) {
	  free_events_.insert(k);
	}
      }
      if (ep.second >= sim_->current_time()) {
	Event_key k= sim_->new_event(ep.second, CGAL_AOS3_TYPENAME EP::U3_event(j_, a,b,c));
	if (k != Event_key() && k != sim_->null_event()) {
	  free_events_.insert(k);
	}
      }
    }
  }
}

/*CGAL_AOS3_TEMPLATE
void Event_visitor CGAL_AOS3_TARG::process_triple(Halfedge_handle h) {
  CGAL_LOG(Log::LOTS, "Processing triple for " << h->opposite()->vertex()->point() 
	   << "--" << h->curve() << "--" << h->vertex()->point() << std::endl);
  if (h->curve().is_arc()
      && h->next()->curve().is_arc()
      && h->next()->next()->curve().is_arc()) {
    Sphere_3_key a= h->curve().key();
    Sphere_3_key b= h->next()->curve().key();
    Sphere_3_key c= h->next()->next()->curve().key();
    
    if (a !=b && b!= c && a != c 
	&& a.is_input() && b.is_input() && c.is_input() 
	&& triples_.find(Sphere_key_utriple(a,b,c)) == triples_.end()) {
      triples_.insert(Sphere_key_utriple(a,b,c));
      Event_pair ep= tr_.intersection_3_events(a,b,c);  
      if (ep.first.is_valid()) {
	CGAL_assertion(ep.first <= ep.second);
	if (a>b) std::swap(a,b);
	if (b>c) std::swap(b,c);
	if (a>b) std::swap(a,b);
     
	if (ep.first >= sim_->current_time()) {
	  Event_key k= sim_->new_event(ep.first, CGAL_AOS3_TYPENAME EP::I3_event(j_, a,b,c));
	  if (k != Event_key() && k != sim_->null_event()) {
	    free_events_.insert(k);
	  }
	}
	if (ep.second >= sim_->current_time()) {
	  Event_key k= sim_->new_event(ep.second, CGAL_AOS3_TYPENAME EP::U3_event(j_, a,b,c));
	  if (k != Event_key() && k != sim_->null_event()) {
	    free_events_.insert(k);
	  }
	}
      }
    }
  }
  }*/

CGAL_AOS3_END_INTERNAL_NAMESPACE
