#include <CGAL/Arrangement_of_spheres_3/Slice.h>
#include <CGAL/Kinetic/Event_base.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_impl.h>
#include <CGAL/IO/Qt_examiner_viewer_2.h>

//typedef CGAL::Kinetic::Event_base<Slice*> Event_base;

struct Slice::Event_base: CGAL::Kinetic::Event_base<Slice*> {
  typedef CGAL::Kinetic::Event_base<Slice*> P;
  Event_base(Slice *sw): P(sw){}
  CGAL::Comparison_result compare_concurrent(Slice::Simulator::Event_key a,
					     Slice::Simulator::Event_key b) const {
    return P::kds()->compare_concurrent(a,b);
  }
};

struct Slice::One_sphere_event: public Event_base {
  typedef Event_base P;
  
  One_sphere_event(Slice *sw, T::Key k): P(sw), k_(k){}

  void process() {
    P::kds()->process_one_sphere_event(k_);
  }
  
  void write(std::ostream &out) const {
    out << k_;
  }

  void audit(Slice::Event_key k) const {
    P::kds()->audit_one_sphere_event(k_, k);
  }

  Slice::T::Key k_;
};

struct Slice::Two_sphere_event: public Event_base {
  typedef Event_base P;
  
  Two_sphere_event(Slice *sw, T::Key a, T::Key b, bool first): P(sw), a_(a), b_(b), first_(first){}

  void process() {
    P::kds()->process_two_sphere_event(a_, b_, first_);
  }
  
  void write(std::ostream &out) const {
    out << a_ << " " << b_ << " (" << first_ << ")"; 
  }

  void audit(Slice::Event_key k) const {
    P::kds()->audit_two_sphere_event(a_, b_, first_, k);
  }

  T::Key a_, b_;
  bool first_;
};


struct Slice::Three_sphere_event: public Event_base {
  typedef Event_base P;
  
  Three_sphere_event(Slice *sw, T::Key a, T::Key b, 
		     T::Key c, bool first): P(sw), a_(a), b_(b), c_(c), first_(first){}

  void process() {
    P::kds()->process_three_sphere_event(a_, b_, c_, first_);
  }
  
  void write(std::ostream &out) const {
    out << a_ << " " << b_ << " " << c_ << " (" << first_ << ")"; 
  }

  void audit(Slice::Event_key k) const {
    P::kds()->audit_three_sphere_event(a_, b_, c_, first_, k);
  }

  T::Key a_, b_, c_;
  bool first_;
};


struct Slice::Rule_event: public Event_base {
  typedef Event_base P;
  
  Rule_event(Slice *sw, T::Key a, 
	     int rule, T::Key o): P(sw), a_(a), ok_(o), rule_(rule){}

  void process() {
    P::kds()->process_rule_collapse_event(a_, rule_, ok_);
  }
  
  void write(std::ostream &out) const {
    out << a_ << " " << rule_ << " " << ok_; 
  }

  void audit(Slice::Event_key k) const {
    P::kds()->audit_rule_collapse_event(a_,  rule_, ok_, k);
  }

  T::Key a_, ok_;
  int rule_;
};


struct Slice::Edge_event: public Event_base {
  typedef Event_base P;
  
  Edge_event(Slice *sw, Halfedge_handle a): P(sw), h_(a){}

  void process() {
    P::kds()->process_edge_collapse_event(h_);
  }
  
  void write(std::ostream &out) const {
    out << h_->vertex()->point() << "--" << h_->curve() << "--"
	<< h_->opposite()->vertex()->point(); 
  }

  void audit(Slice::Event_key k) const {
    P::kds()->audit_edge_collapse_event(h_, k);
  }

  Halfedge_handle h_;
};


Slice::Slice(T tr): t_(tr), sds_(t_.number_of_spheres()),
		    sim_(new Simulator(Simulator::Time(tr.bbox_3().xmin()-1), 
				       Simulator::Time(tr.bbox_3().xmax()+1))),
		    siml_(sim_, this), guil_(NULL){
  /*double d0(tr.bbox_3().xmin()-1); 
  double d1(tr.bbox_3().xmax()+1);
  Simulator::Time t0(tr.bbox_3().xmin()-1); 
  Simulator::Time t1(tr.bbox_3().xmax()+1);*/
  //std::cout << d0 << " " << d1 << " " << t0 << " " << t1 << std::endl;
  //std::cout << *sim_ << std::endl;
  
  double d0(tr.bbox_3().xmin()-1); 
  double d1(tr.bbox_3().xmax()+1);
  Simulator::Time t0(d0); 
  Simulator::Time t1(d1);
  std::cout << d0 << " " << d1 << " " << t0 << " " << t1 << std::endl;
  sim_->set_interval(t0, t1);
  initialize_certificates();
}


void Slice::initialize_certificates() {
  //boost::array<Vertex_handle, 4> dv={{0}};
  //rule_events_.resize(t_.number_of_spheres());
  //CGAL_assertion(rule_certificates_.size() == halfedges_.size());
  for (T::Sphere_key_iterator it = t_.sphere_keys_begin();
       it != t_.sphere_keys_end(); ++it) {
    T::Event_pair ep= t_.sphere_events(*it);
    CGAL_assertion(ep.first <= ep.second);
    if (ep.first >= sim_->current_time()) {
      sim_->new_event(ep.first, One_sphere_event(this, *it));
      //std::cout << *sim_;
    } 
    if (ep.second >= sim_->current_time()) {
      sim_->new_event(ep.second, One_sphere_event(this, *it));
      //std::cout << *sim_;
    }
  }
}


void Slice::audit_one_sphere_event(T::Key k, Event_key ek){
  
}
void Slice::audit_two_sphere_event(T::Key k, T::Key l, bool first, Event_key ek){
  Intersection_2 i2(k,l);
  CGAL_assertion(intersections_2_.find(i2) != intersections_2_.end());
  if (first) {
    CGAL_assertion(intersections_2_[i2].first ==ek);
  } else {
    CGAL_assertion(intersections_2_[i2].second ==ek);
  }
}
void Slice::audit_three_sphere_event(T::Key k, T::Key l, T::Key m, bool first, 
				     Event_key ek){
  Intersection_3 i3(k,l, m);
  CGAL_assertion(intersections_3_.find(i3) != intersections_3_.end());
  if (first) {
    CGAL_assertion(intersections_3_[i3].first ==ek);
  } else {
    CGAL_assertion(intersections_3_[i3].second ==ek);
  }
}
void Slice::audit_rule_collapse_event(T::Key k, int r, T::Key o, Event_key ek){
  Rule_event_rep rep(k,r,o);
  CGAL_assertion(rule_events_.find(rep) != rule_events_.end());
  CGAL_assertion(rule_events_[rep].first == ek 
		 || rule_events_[rep].second==ek);
}
void Slice::audit_edge_collapse_event(Halfedge_handle h, Event_key ek){
  CGAL_assertion(h->event() == ek);
}

void Slice::set_gui(Qt_gui::Handle qt){
  CGAL_precondition(guil_==NULL);
  guil_= new Gui_listener(qt, this);
}

void Slice::process_one_sphere_event(T::Key k) {
  std::cout << "Processing event for sphere " << k 
	    << " at time " << sim_->current_time() << std::endl;
  if (sds_.a_halfedge(k) ==  Halfedge_handle()) {
    insert_sphere(sim_->current_time(), k);
  } else {
    erase_sphere(sim_->current_time(), k);
  }
}

void Slice::process_two_sphere_event(T::Key k, T::Key l, bool f) {
  std::cout << "Two sphere kds " << k << " " << l << " " << f << std::endl;
  Intersection_2 i2(k,l);
  if (!f) {
    intersections_2_.erase(i2);
  } else {
    intersections_2_[i2].first=Event_key();
  }
  // check if centers align
  CGAL::Comparison_result c[2];
  c[0]= t_.compare_sphere_centers_c(k,l, plane_coordinate(0));
  c[1] t_.compare_sphere_centers_c(k,l, plane_coordinate(1));
  if (c[0]== CGAL::EQUAL && c[1] == CGAL::EQUAL) {
    if (t_.compare_sphere_centers_c(k,l, sweep_coordinate()) == CGAL::LARGER) {
      std::swap(k,l);
    }
    sds_.exchange_spheres(k,l);
    sim_->delete_event(intersections_2_[i2].second);
    intersections_2_.erase(i2);
    
    for (unsigned int i=0; i< 4; ++i) {
      Halfedge_handle rh= sds_.rule_halfedge(k, i);
      // now outside
      clean_edge(rh);
      check_edge_collapse(rh);
    }
    for (unsigned int i=0; i< 4; ++i) {
      Halfedge_handle rh= sds_.rule_halfedge(l, i);
      // now outside
      clean_edge(rh);
      if (!(rh->vertex()->point().is_sphere_rule() 
	    && rh->vertex()->point().sphere_key() == k)) {
	check_edge_collapse(rh);
      }
    }
  } else if (c[0]== CGAL::EQUAL || c[1]== CGAL::EQUAL) {
    // figure out which extremums
    const int equal_dir= (c[0]== CGAL::EQUAL)? 0:1;
      const int other_dir= 1-equal_dir;
      // make sure that k is smaller
      if (c[other_dir]== CGAL::LARGER) {
	std::swap(k,l);
      }
    if (f) {
      // has to be outside
    
      int rk= other_dir, rl= other_dir+2;;
      {
	Halfedge_handle rule_k = sds_.rule_halfedge(k, rk);
	Halfedge_handle rule_l = sds_.rule_halfedge(l, rl);
	CGAL_assertion(rule_k->opposite() == rule_l); //otherwise degeneracy
	CGAL_assertion(rule_k->event() == Event_key());
	
	sds_.remove_rule(rule_k);
      }
      Halfedge_handle extr_k = sds_.extremum_halfedge(k, rk);
      Halfedge_handle extr_l = sds_.extremum_halfedge(l, rl);

      // name them for who connects later
      Halfedge_handle rule_vertex_l;
      if (extr_k->next()->key() != k) {
	clean_edge(extr_k->next());
	std::pair<Halfedge_handle, Halfedge_handle> hr
	  = sds_.remove_rule(extra_k->next());
	rule_vertex_l= hr.first;
	CGAL_assertion(hr.second == extr_k);
      }
      Halfedge_handle rule_vertex_k;
      if (extr_l->next()->key() != l) {
	clean_edge(extr_l->next());
	std::pair<Halfedge_handle, Halfedge_handle> hr
	  = sds_.remove_rule(extra_l->next());
	rule_vertex_k= hr.first;
	CGAL_assertion(hr.second == extr_l);
      }
      // maybe replace by two pinch operations
      std::pair<Halfedge_handle, Halfedge_handle> hp= 
	intersect_2(extr_l->next(), extr_l->opposite(),
		    extr_k->next(), extr_k->opposite());
    
      CGAL_assertion(hp.first->curve().key() == k);
      CGAL_assertion(hp.second->curve().key() == l);
      if (rule_vertex_k== Halfedge_handle()) {
	rule_vertex_k= find_rule_vertex(sim_->current_time(),
					hp.first->opposite()->face(),
					Sds::Curve::make_rule(k, rk));
      } 
      sds_.split_face(hp.first->prev(), rule_vertex_k,
		      Sds::Curve::make_rule(k, rk));
      if (rule_vertex_l== Halfedge_handle()) {
	rule_vertex_l= find_rule_vertex(sim_->current_time(),
					hp.second->opposite()->face(),
					Sds::Curve::make_rule(l, rl));
      } 
      sds_.split_face(hp.second->prev(), rule_vertex_l,
		      Sds::Curve::make_rule(l, rl));
      
      // set halfedges for spheres
      sds_.set_extremum_edge(k, hp.first->prev());
      sds_.set_extremum_edge(l, hp.second->prev());
    } else {
      // first figure out which case we are in
      const T::Event_point_3 ep& = sim_->current_time();
      CGAL::Comparison_result c[2];
      c[0]= t_.compare_sphere_center_c(k, ep, T::Coordinate_index(other_dir));
      c[1]= t_.compare_sphere_center_c(l, ep, T::Coordinate_index(other_dir));
      if (c[0] != c[1]) {
	
      } else {
	// yeah!!
	
      }
      // then handle hard case
      // handle easy case by accident
    }

  } else if (f) {
    typedef std::pair<Halfedge_handle, Halfedge_handle> EP;
    EP ep;
    {
      std::vector<Face_handle> shared_faces;
      Halfedge_handle ke, le;
      Halfedge_handle kh= sds_.a_halfedge(k), kh_end=kh;
      do {
	Halfedge_handle oh= kh->opposite()->next()->next();
	do {
	  if (oh->curve().key() == l && oh->curve().is_arc()) {
	    shared_faces.push_back(oh->face());
	    ke= kh->opposite();
	    le= oh;
	  }
	  oh= oh->next();
	} while (oh != kh->opposite());
	
	kh= sds_.next_edge_on_curve(kh);
	CGAL_assertion(kh != Halfedge_handle());
      } while (kh != kh_end);
      // find edges
      CGAL_assertion(!shared_faces.empty());
      if (shared_faces.size() != 1) {
	shared_faces.erase(std::uniq(shared_faces.begin(), shared_faces.end()),
			   shared_faces.end());
	Face_handle fh= locate_point(shared_faces.begin(), shared_faces.end(),
				     ep, key);
      } 
      ep= shared_faces[0];
    }
  } else {
    Face_handle f;
    Halfedge_handle kh= sds_.a_halfedge(k);
    Halfedge_handle khs=kh;
    do {
      if (kh->next()->curve().key() == l 
	  && kh->next()->next() == kh){
	f= kh->face();
	break;
      } else if (kh->opposite()->next()->curve().key() ==l
	  && kh->opposite()->next()->next() == kh->opposite()) {
	f= kh->opposite()->face();
	break;
      }
      kh= sds_.next_edge_on_curve(kh);
      CGAL_assertion(kh != Halfedge_handle());
    } while (kh != khs);
    // find edges
    //CGAL_assertion(!shared_faces.empty());
    unintersect_spheres(sim_->current_time(), f);
  }
}

void Slice::process_three_sphere_event(T::Key k, T::Key l, T::Key m, bool f) {
  std::cout << "Three sphere kds " << k << " " << l << " " << m << " " 
	    << f << std::endl;
  Intersection_3 i3(k,l,m);
  if (!f) {
    intersections_3_.erase(i3);
  } else {
    intersections_3_[i3].first=Event_key();
  }
}

//must maintain hafledges_

void Slice::process_rule_collapse_event(T::Key k, int rule, T::Key o) {
  std::cout << "Rule kds " << k << " " << rule << " " << o << std::endl;
  Rule_event_rep rep(k, rule, o);
  CGAL_precondition(rule_events_.find(rep) != rule_events_.end());
  if (rule_events_[rep].first == Event_key()) {
    rule_events_.erase(rep);
  } else {
    rule_events_[rep].first= Event_key();
  }

  // find halfedge
  Halfedge_handle h= sds_.rule_halfedge(k, rule);
  Halfedge_handle base;
  if (h->next()->next()->next() == h) {
    base= h->next();
  } else {
    if (h->opposite()->next()->next()->next() == h->opposite()) {
      base= h->opposite()->prev();
    }
  }

  if (base != Halfedge_handle()) {
    // we have a uncollapse
    CGAL_assertion(0);
    //uncollapse_edge(sim_->current_time(), h, base);
  } else {
    if (h->next()->next()->curve().key() == o) {
      CGAL_assertion(h->next()->next()->curve().is_arc());
      base= h->next();
    } else {
      base= h->opposite()->prev();
    }
    CGAL_assertion(0);
    
    //collapse_edge(h, base);
  }
  //CGAL_assertion(oh->vertex()->point().is_extremum());
  //CGAL_assertion(oh->curve().key() == k);
  
  
  // see which way it goes
  // erase rule 
  // insert new rule
  // clear event
}

void Slice::process_edge_collapse_event(Halfedge_handle h) {
  std::cout << "Rule edge "; sds_.write(h, std::cout);
  std::cout << std::endl;
  sds_.unset_event(h);

  if (h->curve().is_rule()) {
    Halfedge_handle nh= sds_.next_edge_on_curve(h);
    if (nh == Halfedge_handle()) {
      h=h->opposite();
      nh= sds_.next_edge_on_curve(h);
    }
    if (nh == Halfedge_handle()) {
      Halfedge_handle remove= h->next()->opposite();
      if (remove->curve().is_outward_rule()) {
	remove= h->opposite()->prev();
      }
      CGAL_precondition(!remove->curve().is_outward_rule());
      CGAL_precondition(remove->curve().key() != h->curve().key());
      /*Halfedge_handle nr=*/ rotate_rule(sim_->current_time(), remove);
    }
    
  }
  
  Halfedge_handle p= sds_.cross_edge(h)->opposite();
  CGAL_assertion(p != Halfedge_handle());
  if (!p->curve().is_rule()) {
    // h=h->opposite();
    p= sds_.cross_edge(h->opposite())->opposite();
  }
  // simplify?
  if (h->opposite()->face() == p->face() 
      || h->opposite()->face() == p->opposite()->face()) {
    h=h->opposite();
  }

  CGAL_assertion(p->curve().is_rule());
  //CGAL_assertion(p->vertex() == h->vertex());

  clean_edge(h->prev());
  clean_edge(h->next());
  clean_edge(h->opposite()->prev());
  clean_edge(h->opposite()->next());
  bool left=false;
  if (p==h->prev()) {
    left=true;
    clean_edge(p->opposite()->prev());
  } else {
    clean_edge(p->next());
  }
  
  //Halfedge_handle rule, prev, next;
  
  // rule->vertex() is constant
  // all faces are constant

  sds_.exchange_vertices(h, p);
  
  
  Qt_examiner_viewer_2 *qt= new Qt_examiner_viewer_2();
  draw_rz(qt, CGAL::to_double(sim_->current_time()) + .1);
  qt->show_everything();
  qt->show();
  

  check_edge_face(h);
  check_edge_collapse(h);
  check_edge_collapse(h->prev());
  check_edge_collapse(h->next());
  if (left) {
    check_edge_collapse(h->next()->opposite()->next());
    check_reduced_face(h->prev()->opposite()->face());
  } else {
    check_edge_collapse(h->prev()->opposite()->prev());
    check_reduced_face(h->next()->opposite()->face());
  }
  //audit();
}



void Slice::check_edge_pair(Halfedge_handle h,
			    Halfedge_handle k) {
  if (h->curve().key() == k->curve().key() 
      || !h->curve().is_arc() || !k->curve().is_arc()
      || h->curve().quadrant() & k->curve().quadrant() != 0 ) return;

  Intersection_2 i2(h->curve().key(), k->curve().key());
  if (intersections_2_.find(i2) == intersections_2_.end()) {
    T::Event_pair ep= t_.intersection_2_events(i2.first(), i2.second());
    if (ep.first.is_valid()) {
      Event_key ek0= sim_->new_event(ep.first, Two_sphere_event(this,
								i2.first(),
								i2.second(),
								true));
      Event_key ek1= sim_->new_event(ep.second, Two_sphere_event(this,
								 i2.first(),
								 i2.second(),
								 false));
      //std::cout << *sim_;
      intersections_2_[i2]= Event_key_pair(ek0, ek1);
    } else {
      intersections_2_[i2]= Event_key_pair();
    }
  }
}




void Slice::check_edge_face(Halfedge_handle h) {
  Face_handle f= h->face();
  Halfedge_handle c= f->halfedge();
  do {
    check_edge_pair(h, c);
    c= c->next();
  } while (c != f->halfedge());
}




void Slice::check_merged_faces(Face_handle f, Face_handle g) {
  Halfedge_handle cf= f->halfedge();
  do {
    Halfedge_handle cg= g->halfedge();
    do {
      check_edge_pair(cg, cf);
      cg= cg->next();
    } while (cg != g->halfedge());
    
    cf= cf->next();
  } while (cf != f->halfedge());
}





void Slice::check_edge_collapse(Halfedge_handle h) {
  // complicated a bit
  // elimitate two SS endpoints
  std::cout << "Checking edge collapse for " << h->vertex()->point() 
	    << "--" << h->curve() << "--" << h->opposite()->vertex()->point()
	    << std::endl;
  CGAL_precondition(sds_.event(h)== Event_key());
  if (!h->curve().is_finite()) return;
  // an arc connecting spheres will be checked elsewhere
  if (h->vertex()->point().is_sphere_sphere() 
      && h->opposite()->vertex()->point().is_sphere_sphere()) return;
  // an arc on a sphere with its two rules (but they might be called
  // something else
  if (h->vertex()->point().is_sphere_extremum() 
      && h->opposite()->vertex()->point().is_sphere_extremum()) return;
  // rule collapse will be detected separately
  if (h->curve().is_rule() && h->vertex()->point().is_sphere_rule()
      && h->opposite()->vertex()->point().is_sphere_rule()){
    if (h->vertex()->point().is_sphere_extremum() 
	&& h->opposite()->vertex()->point().is_sphere_extremum()) {
      // degenerate edge check separately, ick.
      CGAL_assertion(0);
    }
    return;
  }

  Halfedge_handle rh= sds_.cross_edge(h)->opposite();
  Halfedge_handle oh= sds_.cross_edge(h->opposite())->opposite();
  if (!rh->curve().is_rule()) std::swap(rh, oh);
  
  if (!oh->curve().is_finite()) return;
  if (!rh->curve().is_finite()) return;
  std::cout << "Found " << rh->curve() 
	    << ", " << h->curve() << "--" << oh->curve()
	    << std::endl;
  CGAL_assertion(h->curve() != rh->curve());
  CGAL_assertion(oh->curve() != rh->curve());
  CGAL_assertion(h->curve() != oh->curve());
  CGAL_assertion(rh->curve().is_rule());
  if (h->curve().key() == rh->curve().key()
      && h->curve().key() == oh->curve().key()) return;
  if (h->curve().key() == rh->curve().key() && !oh->curve().is_rule()) {
    // extremum
    Rule_event_rep re(rh->curve().key(), rh->curve().rule_index(), oh->curve().key());
    if (rule_events_.find(re) == rule_events_.end()) {
      T::Event_pair ep
	= t_.sphere_intersect_extremum_events(re.rule_key(), 
					      re.rule_coordinate(),
					      re.other_key());
      if (ep.first.is_valid() && ep.first != ep.second) {
	// not handling degeneracies
	Event_key ek0, ek1;
	if (ep.first > sim_->current_time()) {
	  ek0= sim_->new_event(ep.first, Rule_event(this,
						    re.rule_key(), 
						    re.rule_index(),
						    re.other_key()));
	} 
	if (ep.second > sim_->current_time()) {
	  ek1= sim_->new_event(ep.second, Rule_event(this,
						     re.rule_key(), 
						     re.rule_index(),
						     re.other_key()));
	}
	rule_events_[re]= Event_key_pair(ek0, ek1);
	//std::cout << *sim_;
	std::cout << "Created rule events " << ek0 << " " << ek1 << std::endl;
      } else {
	rule_events_[re]= Event_key_pair();
      } 
    }
  } else if (h->curve().is_rule()) {
    //SR cross RR
    T::Key ks[2];
    ks[project(h->curve().constant_coordinate())]= h->curve().key();
    ks[project(rh->curve().constant_coordinate())]= rh->curve().key();
    T::Event_pair ep= t_.sphere_intersect_rule_rule_events(oh->curve().key(),
							   ks[0],
							   ks[1]);
    Event_key ek= sim_->null_event();
    // handle degeneracy here by looking at whether I am inside already
    if (ep.first.is_valid() && ep.first != ep.second) {
      if (ep.first  > sim_->current_time()) {
	ek = sim_->new_event(ep.first, Edge_event(this, h));
      } else if (ep.second > sim_->current_time()) {
	ek= sim_->new_event(ep.second, Edge_event(this,h));
      }
    }
    // std::cout << *sim_;
    sds_.set_event(h, ek);
    std::cout << "Created SR RR collapse event " << ek << std::endl;
  } else if (oh->curve().is_rule()) {
    // parallel rules
    if (oh->curve().constant_coordinate() == rh->curve().constant_coordinate()) {
      return;
    }
    // SR cross SR
    T::Key ks[2];
    ks[project(oh->curve().constant_coordinate())]= oh->curve().key();
    ks[project(rh->curve().constant_coordinate())]= rh->curve().key();
    T::Event_pair ep= t_.sphere_intersect_rule_rule_events(h->curve().key(),
							   ks[0],
							   ks[1]);
    Event_key ek= sim_->null_event();
    // handle degeneracy here by looking at whether I am inside already
    if (ep.first.is_valid()  && ep.first != ep.second) {
      if (ep.first  > sim_->current_time()) {
	ek = sim_->new_event(ep.first, Edge_event(this, h));
      } else if (ep.second > sim_->current_time()) {
	ek= sim_->new_event(ep.second, Edge_event(this,h));
      }
    }
    //std::cout << *sim_;
    sds_.set_event(h, ek);
    std::cout << "Created SR SR collapse event " << ek << std::endl;
  } else {
    // SS cross SR
    T::Event_pair ep= t_.circle_cross_rule_events(h->curve().key(),
						  oh->curve().key(),
						  rh->curve().key(),
						  rh->curve().constant_coordinate());
    Event_key ek= sim_->null_event();
    // handle degeneracy here by looking at whether I am inside already
    if (ep.first.is_valid() && ep.first != ep.second) {
      if (ep.first  > sim_->current_time()) {
	ek = sim_->new_event(ep.first, Edge_event(this, h));
      } else if (ep.second > sim_->current_time()) {
	ek= sim_->new_event(ep.second, Edge_event(this,h));
      }
    }
    //std::cout << *sim_;
    sds_.set_event(h, ek);
    std::cout << "Created SS SR collapse event " << ek << std::endl;
  }
}









void Slice::check_reduced_face(Face_handle f) {
  if (f->halfedge()->next()->next()->next() == f->halfedge()) {
    Halfedge_handle r= f->halfedge();
    do {
      if (r->curve().is_rule()) break;
      r=r->next();
    } while (r != f->halfedge());

    if (r->curve().is_rule()) {
      T::Key rk= r->curve().key();
      int ri= r->curve().rule_index();
      T::Key ok;
      if (r->next()->curve().key() != rk) {
	ok = r->next()->curve().key();
      } else {
	ok = r->prev()->curve().key();
      }
      CGAL_assertion(ok != rk);
      
      Rule_event_rep re(rk, ri, ok);
      if (rule_events_.find(re) == rule_events_.end()) {
	T::Event_pair ep
	  = t_.sphere_intersect_extremum_events(re.rule_key(), 
						re.rule_coordinate(),
						re.other_key());
	if (ep.first.is_valid()) {
	  Event_key ek0= sim_->new_event(ep.first, Rule_event(this,
							      rk, ri, ok));
	  Event_key ek1= sim_->new_event(ep.second, Rule_event(this,
							       rk, ri, ok));
	  rule_events_[re]= Event_key_pair(ek0, ek1);
	  std::cout << *sim_;
	} else {
	  rule_events_[re]= Event_key_pair();
	}
      }
    } else if (f->halfedge()->curve().key() !=
	      f->halfedge()->next()->curve().key()
	      && f->halfedge()->next()->next()->curve().key() 
	       != f->halfedge()->curve().key()
	      && f->halfedge()->next()->curve().key()
	       != f->halfedge()->next()->next()->curve().key()) {
      Intersection_3 i3(f->halfedge()->curve().key(),
			f->halfedge()->next()->curve().key(),
			f->halfedge()->next()->next()->curve().key());
      if (intersections_3_.find(i3) == intersections_3_.end()) {
	T::Event_pair ep= t_.intersection_3_events(i3.first(), i3.second(), i3.third());
	if (ep.first.is_valid()) {
	  
	  Event_key ek0= sim_->new_event(ep.first, Three_sphere_event(this,
								      i3.first(),
								      i3.second(),
								      i3.third(),
								      true));
	  Event_key ek1= sim_->new_event(ep.second, Three_sphere_event(this,
								       i3.first(),
								       i3.second(),
								       i3.third(),
								       false));
	  std::cout << *sim_;
	  intersections_3_[i3]= Event_key_pair(ek0, ek1);
	} else {
	  intersections_3_[i3]= Event_key_pair();
	}
      } 
    }
  }
}




void Slice::clean_edge(Halfedge_handle h) {
  if (h->event() != Event_key()) {
    sim_->delete_event(h->event());
    sds_.unset_event(h);
  }
}

/*
  void Slice::cross_vertex_event(Halfedge_handle c) {
  
  // degeneracies will be a pain
  CGAL_assertion(sds_.degree(c->opposite()->vertex())==3); 
  clean_edge(c);
  clean_edge(c->prev());
  clean_edge(c->next());
  clean_edge(c->prev()->opposite()->prev());
  
  Halfedge_handle rule, prev, next;
  
  // rule->vertex() is constant
  // all faces are constant

  sds_.exchange_vertices(c);
  

  on_new_edge(c->prev());
  on_new_edge(c->next());
  on_new_edge(c);
  on_new_edge(c->next()->opposite()->next());
  }*/




void Slice::audit_events() const {
  /*std::cout << "Auditing certificates..." << std::flush;
    boost::array<bool, 4> da={{false,false,false,false}};
    std::vector<boost::array<bool, 4> > rule_certs(t_.number_of_spheres(), 
    da);
    for (Sds::Halfedge_const_iterator vit= sds_.halfedges_begin();
    vit != sds_.halfedges_end(); ++vit) {
    if (vit->curve().is_rule() && vit->curve().is_finite()) {
    Halfedge_const_handle h= vit;
    do {
    Sds::Point vp=h->vertex()->point(); 
    Sds::Point vop= h->opposite()->vertex()->point();
    Sds::Curve vc= h->curve();
    if (vp.is_sphere_rule()
    && h->curve().key() == h->next()->curve().key()) {
    //Sds::Curve oc= vop.other_curve(vc);
    T::Key k=vc.key();
    int ii= k.input_index();
    //int sz= rule_certificates_.size();
    Vertex_const_handle tvh= rule_events_[ii][vc.rule_index()];
    CGAL_assertion(tvh == h->opposite()->vertex());
    rule_certs[vc.key().input_index()][vc.rule_index()]=true;
    }
    h= h->opposite();
    } while (h != vit);

    }
    }

    for (unsigned int i=0; i< rule_certs.size(); ++i){
    for (unsigned int j=0; j< 4; ++j){
    if (!rule_certs[i][j]) {
    CGAL_assertion(rule_certificates_[i][j] == Vertex_handle());
    }
    }
    }*/
  /*for (Sds::Halfedge_const_iterator h= sds_.halfedges_begin();
    h != sds_.halfedges_end(); ++h) {
    if (h->curve().is_rule()) {
    // if both SR then they may or may not have a certi
    // if one RR and one SR, check if they intersect if so check for certificate
    if (h->vertex()->point().is_rule_rule()) {
    h=h->opposite();
    }
    if (h->vertex()->point().is_sphere_rule() && h->opposite()->vertex()->point().is_rule_rule()) {
    // check if they intersect, if so, must have event
    Sds::Curve r= h->opposite()->vertex()->point().rule(1-h->vertex()->point().rule_coordinate());
    if (t_.sphere_intersects_rule(h->vertex()->point().sphere_key(),
    r.key(),
    r.constant_coordinate())) {
    T::Event_pair ep= t_.intersection_events(h->vertex()->point().sphere_key(),
    r.key(),
    r.constant_coordinate());
    if (ep.second >= sim_->current_time()) {
    CGAL_assertion(h->event() != Simulator::Event_key());
    } else {
    CGAL_assertion(h->event() == Simulator::Event_key());
    }
    }
    }
    } else {
    if (!h->vertex()->point().is_sphere_rule()) {
    h=h->opposite();
    }
      
    }
    if (intersects(h->vertex()->point(), h->opposite()->vertex()->point())) {
    T::Event_pair ep= intersection_events(h->vertex()->point(), h->opposite()->vertex()->point());
    if (ep.second >= sim_->current_time()){
    CGAL_assertion(h->event() != Simulator::Event_key());
    } else {
    CGAL_assertion(h->event() == Simulator::Event_key());
    }
    }
    }*/


  /*for (Sds::Face_const_iterator vit= sds_.faces_begin();
    vit != sds_.faces_end(); ++vit) {
    Sds::Halfedge_const_handle h= vit->halfedge();
    do {
    Sds::Halfedge_const_handle h2= vit->halfedge();
    do {
    if (h2->curve().is_finite() && h->curve().is_finite() 
    && h2->curve().is_arc() && h->curve().is_arc() 
    && h2->curve().key() != h->curve().key() 
    && !(h2->curve().quadrant() & h->curve().quadrant())
    && t_.intersects(h2->curve().key(), h->curve().key())) {
    Intersection_2 i2(h->curve(),
    h2->curve());
    CGAL_assertion(certificates_2_.find(i2) 
    != certificates_2_.end());
    //int num=certificates_2_.find(i2)->second;
    //CGAL_assertion(num<=2 && num >=0);
    }
    h2 = h2->next();
    } while (h2 != h);
    h=h->next();
    } while (h != vit->halfedge());
    }*/
#if 0
  for (Sds::Face_const_iterator vit= sds_.faces_begin();
       vit != sds_.faces_end(); ++vit) {
    if ( vit->halfedge()->next()->vertex() 
	 ==  vit->halfedge()->prev()->opposite()->vertex()) {
      // it has size 3
      Sds::Halfedge_const_handle h= vit->halfedge();
      if (h->curve().is_arc() && h->next()->curve().is_arc()
	  && h->prev()->curve().is_arc()
	  && t_.intersects(h->curve().key(),
			   h->next()->curve().key(),
			   h->prev()->curve().key())) {
	Intersection_3 i3(h->curve(),
			  h->next()->curve(),
			  h->prev()->curve());
	/*if (!h->curve().is_inside() && !h->next()->curve().is_inside()
	  && !h->prev()->curve().is_inside()) {*/
	CGAL_assertion(certificates_3_.find(i3)
		       != certificates_3_.end());
	//CGAL_assertion(certificates_3_.find(i3)->second==2);
	/*} else {
	  CGAL_assertion(certificates_3_.find(i3)
	  != certificates_3_.end());
	  //CGAL_assertion(certificates_3_.find(i3)->second==1);
	  }*/
      }
    }
  }
  std::cout << "done." << std::endl;
#endif
}




void Slice::draw_events_rz(CGAL::Qt_widget *qtv, NT z) {
  //t_.set_temp_sphere(T::Sphere_3(T::Point_3(0,0,z), 0));


  *qtv << CGAL::Color(255, 155, 155);
  for (Intersections_2::const_iterator it= intersections_2_.begin();
       it != intersections_2_.end(); ++it){
    if (it->second.first != Event_key()) {
      Event_key k= it->second.first;
      T::Event_point_3 ep= sim_->event_time(k);
      DT::Point_2 dp(ep.approximate_coordinate(plane_coordinate(0)), 
		     ep.approximate_coordinate(plane_coordinate(1)));
      *qtv << dp;
      
      T::Key a= sim_->event<Two_sphere_event>(k).a_;
      T::Key b= sim_->event<Two_sphere_event>(k).b_;
      
      /*std::ostringstream oss;
	oss << a << " " << b;
	*qtv << oss.str().c_str();*/
    }
    if (it->second.second != Event_key()) {
      Event_key k= it->second.second;
      T::Event_point_3 ep= sim_->event_time(k);
      DT::Point_2 dp(ep.approximate_coordinate(plane_coordinate(0)), 
		     ep.approximate_coordinate(plane_coordinate(1)));
      *qtv << dp;
      
      T::Key a= sim_->event<Two_sphere_event>(k).a_;
      T::Key b= sim_->event<Two_sphere_event>(k).b_;
      
      /*std::ostringstream oss;
	oss << a << " " << b;
	*qtv << oss.str();*/
    }
  }

  *qtv << CGAL::Color(255, 155, 155);
  for (Intersections_3::const_iterator it= intersections_3_.begin();
       it != intersections_3_.end(); ++it) {
    if (it->second.first != Event_key()) {
      Event_key k= it->second.first;
      T::Event_point_3 ep= sim_->event_time(k);
      DT::Point_2 dp(ep.approximate_coordinate(plane_coordinate(0)), 
		     ep.approximate_coordinate(plane_coordinate(1)));
      *qtv << dp;
      
      /*T::Key a= sim_->event<Three_sphere_event>(k).a_;
	T::Key b= sim_->event<Three_sphere_event>(k).b_;
	T::Key c= sim_->event<Three_sphere_event>(k).c_;
      
	std::ostringstream oss;
	oss << a << " " << b << " " << c;
	*qtv << oss.str();*/
    }
    if (it->second.second != Event_key()){
      Event_key k= it->second.second;
      T::Event_point_3 ep= sim_->event_time(k);
      DT::Point_2 dp(ep.approximate_coordinate(plane_coordinate(0)), 
		     ep.approximate_coordinate(plane_coordinate(1)));
      *qtv << dp;
      
      /*T::Key a= sim_->event<Three_sphere_event>(k).a_;
	T::Key b= sim_->event<Three_sphere_event>(k).b_;
	T::Key c= sim_->event<Three_sphere_event>(k).c_;
      
	std::ostringstream oss;
	oss << a << " " << b << " " << c;
	*qtv << oss.str();*/
    }
  }
  //std::set<Intersection_3> certificates_3_;

  
}




CGAL::Comparison_result Slice::compare_concurrent(Simulator::Event_key k0,
			       Simulator::Event_key k1) const {
  const T::Event_point_3& p0= sim_->event_time(k0);
  const T::Event_point_3& p1= sim_->event_time(k1);
  CGAL::Comparison_result ret= p0.compare(p1, plane_coordinate(0));
  if (ret != CGAL::EQUAL) {
    return ret; 
  } else {
    return p0.compare(p1, plane_coordinate(1));
  }
}
