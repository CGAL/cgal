#include <CGAL/Arrangement_of_spheres_3/Slice.h>
#include <CGAL/Kinetic/Event_base.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_impl.h>
#include <CGAL/IO/Qt_examiner_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Slice_impl.h>
#include <boost/iterator/indirect_iterator.hpp>

//typedef CGAL::Kinetic::Event_base<Slice*> Event_base;

struct Slice::Event_base: CGAL::Kinetic::Event_base<Slice*> {
  //enum Type {INSERT=0, DELETE, TWO, THREE, RULE, EDGE}
  typedef CGAL::Kinetic::Event_base<Slice*> P;
  Event_base(Slice *sw/*, Type t*/): P(sw){}
  CGAL::Comparison_result compare_concurrent(Slice::Simulator::Event_key a,
					     Slice::Simulator::Event_key b) const {
    return P::kds()->compare_concurrent(a,b);
  }
  /*Type type() const {
    return type_;
    }*/
  //Type type_;
};

struct Slice::One_sphere_event: public Event_base {
  typedef Event_base P;
  
  One_sphere_event(Slice *sw, T::Key k): P(sw/*,
							    insert:T::INSERT?T::DELETE*/), k_(k){}

  void process() {
    P::kds()->process_one_sphere_event(k_);
  }
  
  void write(std::ostream &out) const {
    out << "I/R " << k_;
  }

  void audit(Slice::Event_key k) const {
    P::kds()->audit_one_sphere_event(k_, k);
  }

  Slice::T::Key k_;
};

struct Slice::Two_sphere_event: public Event_base {
  typedef Event_base P;
  
  Two_sphere_event(Slice *sw, T::Key a, T::Key b, bool first): P(sw/*, T::TWO*/),
							       a_(a), b_(b), first_(first){}

  void process() {
    P::kds()->process_two_sphere_event(a_, b_, first_);
  }
  
  void write(std::ostream &out) const {
    out << "I/U" << a_ << " " << b_ << " (" << first_ << ")"; 
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
		     T::Key c, bool first): P(sw/*, T::THREE*/),
					    a_(a), b_(b), c_(c), first_(first){}

  void process() {
    P::kds()->process_three_sphere_event(a_, b_, c_, first_);
  }
  
  void write(std::ostream &out) const {
    out << "P " << a_ << " " << b_ << " " << c_ << " (" << first_ << ")"; 
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
	     Rule_direction rule, T::Key o): P(sw/*, T::RULE*/), a_(a), ok_(o), rule_(rule){}

  void process() {
    P::kds()->process_rule_collapse_event(a_, rule_, ok_);
  }
  
  void write(std::ostream &out) const {
    out << "Rule " << a_ << " " << rule_ << " " << ok_; 
  }

  void audit(Slice::Event_key k) const {
    P::kds()->audit_rule_collapse_event(a_,  rule_, ok_, k);
  }

  T::Key a_, ok_;
  Rule_direction rule_;
};


struct Slice::Edge_event: public Event_base {
  typedef Event_base P;
  
  Edge_event(Slice *sw, Halfedge_handle a): P(sw/*, T::EDGE*/), h_(a){}

  void process() {
    P::kds()->process_vertex_crossing_event(h_);
  }
  
  void write(std::ostream &out) const {
    out << "Edge " << h_->vertex()->point() << "--" << h_->curve() << "--"
	<< h_->opposite()->vertex()->point(); 
  }

  void audit(Slice::Event_key k) const {
    P::kds()->audit_edge_collapse_event(h_, k);
  }

  Halfedge_handle h_;
};


Slice::Slice(T tr): t_(tr), sds_(t_.number_of_spheres()),
		    sim_(new Simulator(Simulator::Time(tr.bbox_3().xmin()-1))),
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
  //sim_->set_interval(t0, t1);
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
void Slice::audit_rule_collapse_event(T::Key k, Rule_direction r, 
				      T::Key o, Event_key ek){
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

bool Slice::has_degeneracy() const {
  if (sim_->empty()) return false;
  if (sim_->next_event_time().compare(sim_->current_time(),
				      sweep_coordinate()) == CGAL::EQUAL
      && sim_->next_event_time().compare(sim_->current_time(),
					 plane_coordinate(0)) == CGAL::EQUAL
      && sim_->next_event_time().compare(sim_->current_time(),
					 plane_coordinate(0)) == CGAL::EQUAL) {
    std::cerr << "Degeneracy at " << sim_->current_time() << std::endl;
    //process_degeneracy(sim_->current_time(), );
    return true;
  } else {
    return false;
  }
}


void Slice::process_degeneracy() {
  T::Sphere_point_3 ct= sim_->current_time();
  std::cerr << "Degenerate events: \n";
  while (sim_->next_event_time() == ct) {
    std::cerr << sim_->next_event() << std::endl;;
    sim_->delete_event(sim_->next_event());
  }
  std::cerr << std::endl;
  CGAL_assertion(0);
}

void Slice::rebuild_degenerate(const T::Sphere_point_3 &pt,
			       Vertex_handle vh) {
  CGAL_assertion(0);
}


void Slice::process_one_sphere_event(T::Key k) {
  if (has_degeneracy()) {
    sim_->new_event(sim_->current_time(), One_sphere_event(this, k));
    process_degeneracy();
    return;
  }
  std::cout << "Processing event for sphere " << k 
	    << " at time " << sim_->current_time() << std::endl;
  if (sds_.a_halfedge(k) ==  Halfedge_handle()) {
    insert_sphere(sim_->current_time(), k);
  } else {
    erase_sphere(sim_->current_time(), k);
  }
}



void Slice::process_equal_centers_two_sphere_event(T::Key k, T::Key l, bool f) {;
  Intersection_2 i2(k,l);
  if (t_.compare_sphere_centers_c(k,l, sweep_coordinate()) == CGAL::LARGER) {
    std::swap(k,l);
  }
  exchange_circles(k,l);
  sim_->delete_event(intersections_2_[i2].second);
  intersections_2_[i2].first= Event_key();
  intersections_2_[i2].second= Event_key();
  /*{
    Halfedge_handle kc= sds_.a_halfedge(k), ke=kc;
    do {
      kc->curve().set_key(l);
      kc= sds_.next_edge_on_curve(kc);
    } while (kc != ke);
  }
  {
    Halfedge_handle lc= sds_.a_halfedge(l), le=lc;
    do {
      lc->curve().set_key(k);
      lc= sds_.next_edge_on_curve(lc);
    } while (lc != le);
    }*/
  Halfedge_handle keh[4];
  for (unsigned int i=0; i< 4; ++i) {
    Halfedge_handle rh= sds_.rule_halfedge(k, Rule_direction(i));
    keh[i]= rh;
    // now outside
    clean_edge(rh);
    check_edge_collapse(rh);
  }
  for (unsigned int i=0; i< 4; ++i) {
    Halfedge_handle rh= sds_.rule_halfedge(l, Rule_direction(i));
    // now outside
    clean_edge(rh);
    if (!(rh->vertex()->point().is_sphere_rule() 
	  && rh->vertex()->point().sphere_key() == k)) {
      check_edge_collapse(rh);
    }
  }
}



void Slice::process_aligned_centers_two_sphere_event(T::Key k, T::Key l, bool f,
						     CGAL::Comparison_result c[2]) {
  // figure out which extremums
  Intersection_2 i2(k,l);
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
      Halfedge_handle rule_k = sds_.rule_halfedge(k, Rule_direction(rk));
      Halfedge_handle rule_l = sds_.rule_halfedge(l, Rule_direction(rl));
      CGAL_assertion(rule_k->opposite() == rule_l); //otherwise degeneracy
      CGAL_assertion(rule_k->event() == Event_key());
      
      merge_faces(rule_k);
    }
    Halfedge_handle extr_k = sds_.extremum_halfedge(k, Rule_direction(rk));
    Halfedge_handle extr_l = sds_.extremum_halfedge(l, Rule_direction(rl));
      
    // name them for who connects later
    Halfedge_handle rule_vertex_l;
    if (extr_k->next()->curve().key() != k) {
      rule_vertex_l= extr_k->next()->opposite()->prev();;
      merge_faces(extr_k->next());
      
      //CGAL_assertion(hr.second == extr_k);
    }
    Halfedge_handle rule_vertex_k;
    if (extr_l->next()->curve().key() != l) {
      rule_vertex_k= extr_l->next()->opposite()->prev();;
      merge_faces(extr_l->next());
    }
    // maybe replace by two pinch operations
    clean_edge(extr_l->opposite());
    clean_edge(extr_k->opposite());
    CGAL_assertion(0);
    std::pair<Halfedge_handle, Halfedge_handle> hp;
    /*=
      sds_.intersect_2(extr_l->next(), extr_l->opposite(),
      extr_k->next(), extr_k->opposite());*/
    // clean edges and create certificates
    CGAL_postcondition(hp.first->face() == hp.second->face());
    CGAL_postcondition(hp.first->prev()->prev()== hp.second);
    CGAL_postcondition(hp.second->prev()->prev()== hp.first);

    check_edge_collapse(hp.first->opposite()->prev());
    check_edge_collapse(hp.first->opposite()->prev()->opposite()->prev());
    check_edge_collapse(hp.second->opposite()->prev());
    check_edge_collapse(hp.second->opposite()->prev()->opposite()->prev());
    
    CGAL_assertion(hp.first->curve().key() == k);
    CGAL_assertion(hp.second->curve().key() == l);
    if (rule_vertex_k== Halfedge_handle()) {
      rule_vertex_k= find_rule_vertex(sim_->current_time(),
				      hp.first->opposite()->face(),
				      Sds::Curve::make_rule(k, Rule_direction(rk)));
    } 
    sds_.split_face(Sds::Curve::make_rule(k, Rule_direction(rk)),
		    hp.first->prev(), rule_vertex_k);
    if (rule_vertex_l== Halfedge_handle()) {
      rule_vertex_l= find_rule_vertex(sim_->current_time(),
				      hp.second->opposite()->face(),
				      Sds::Curve::make_rule(l, Rule_direction(rl)));
    } 
    sds_.split_face(Sds::Curve::make_rule(l, Rule_direction(rl)),
		    hp.second->prev(), rule_vertex_l);
      
    // set halfedges for spheres
    CGAL_assertion(0);
    //sds_.set_extremum_halfedge(hp.first->prev());
    //sds_.set_extremum_halfedge(hp.second->prev());
  } else {
    // first figure out which case we are in
    const T::Event_point_3 &ep = sim_->current_time();
    CGAL::Comparison_result c[2];
    c[0]= t_.compare_sphere_center_c(k, ep, T::Coordinate_index(other_dir));
    c[1]= t_.compare_sphere_center_c(l, ep, T::Coordinate_index(other_dir));
    CGAL_assertion(0);
    if (c[0] != c[1]) {
	
    } else {
      // yeah!!
	
    }
    // then handle hard case
    // handle easy case by accident
  }
}


void Slice::process_intersect_event(T::Key k, T::Key l) {
  typedef std::pair<Halfedge_handle, Halfedge_handle> EP;
  EP ep;
  Face_handle fh;
  {
    std::vector<Face_handle> shared_faces;
    typedef  std::vector<Face_handle>::iterator SFIT;
    // Halfedge_handle ke, le;
    Halfedge_handle kh= sds_.a_halfedge(k), kh_end=kh;
    do {
      {
	Halfedge_handle oh= kh->opposite()->next()->next();
	do {
	  if (oh->curve().key() == l && oh->curve().is_arc()) {
	    shared_faces.push_back(oh->face());
	    //ke= kh->opposite();
	    //le= oh;
	  }
	  oh= oh->next();
	} while (oh != kh->opposite());
      }
      {
	Halfedge_handle oh= kh->next()->next();
	do {
	  if (oh->curve().key() == l && oh->curve().is_arc()) {
	    shared_faces.push_back(oh->face());
	    //ke= kh->opposite();
	    //le= oh;
	  }
	  oh= oh->next();
	} while (oh != kh);
      }
      kh= sds_.next_edge_on_curve(kh);
      CGAL_assertion(kh != Halfedge_handle());
    } while (kh != kh_end);
    // find edges
    CGAL_assertion(!shared_faces.empty());
    if (shared_faces.size() != 1) {
      std::cout << "Have several shared faces " << std::endl;
      shared_faces.erase(std::unique(shared_faces.begin(), shared_faces.end()),
			 shared_faces.end());
      for (unsigned int i=0; i< shared_faces.size(); ++i) {
	sds_.write_face(shared_faces[i]->halfedge(), std::cout) << std::endl;
      }
      fh= locate_point(SFIT(shared_faces.begin()), SFIT(shared_faces.end()),
		       sim_->current_time());
      std::cout << "Got ";
      sds_.write_face(fh->halfedge(), std::cout) << std::endl;
    } else {
      fh= shared_faces[0];
    }
  }
  Halfedge_handle ha;
  try {
    ha= shoot_rule(sim_->current_time(), fh, sim_->current_time(),
		   Rule_direction(1));
  } catch (On_vertex_exception v) {
    if (v.vertex_handle()->point().is_sphere_rule()) {
      ha=v.vertex_handle()->halfedge();
      while (ha->face() != fh && !ha->curve().is_arc()) {
	if (ha->opposite()->face() == fh) {
	  ha=ha->opposite(); break;
	}
	ha= ha->next()->opposite();
      }
    } else {
      CGAL_assertion(0);
    }
  }
  Halfedge_handle hb;
  try {
    hb= shoot_rule(sim_->current_time(), fh, sim_->current_time(),
		   Rule_direction(3));
  } catch (On_vertex_exception v) {
    if (v.vertex_handle()->point().is_sphere_rule()) {
      hb=v.vertex_handle()->halfedge();
      while (hb->face() != fh && !hb->curve().is_arc()) {
	if (hb->opposite()->face() == fh) {
	  hb=hb->opposite(); break;
	}
	hb= ha->next()->opposite();
      }
    } else {
      CGAL_assertion(0);
    }
  }

  CGAL_assertion(ha->curve().key() != hb->curve().key());
  CGAL_assertion(ha->curve().key() == k
		 || ha->curve().key() == l );
  CGAL_assertion(hb->curve().key() == k
		 || hb->curve().key() == l );
  Face_handle f= intersect_spheres(sim_->current_time(),
				   ha,hb);
}

void Slice::process_unintersect_event(T::Key k, T::Key l) {
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
  unintersect_spheres(sim_->current_time(), f->halfedge(), f->halfedge()->next());
}


void Slice::process_two_sphere_event(T::Key k, T::Key l, bool f) {
  if (has_degeneracy()) {
    sim_->new_event(sim_->current_time(), Two_sphere_event(this, k, l, f));
    process_degeneracy();
    return;
  }
  std::cout << "Two sphere kds " << k << " " << l << " " << f << std::endl;
  Intersection_2 i2(k,l);
  if (!f) {
    //intersections_2_.erase(i2);
    intersections_2_[i2].second=Event_key();
  } else {
    intersections_2_[i2].first=Event_key();
  }
  // check if centers align
  CGAL::Comparison_result c[2];
  c[0]= t_.compare_sphere_centers_c(k,l, plane_coordinate(0));
  c[1]= t_.compare_sphere_centers_c(k,l, plane_coordinate(1));
  if (c[0]== CGAL::EQUAL && c[1] == CGAL::EQUAL) {
    process_equal_centers_two_sphere_event(k,l,f);
  } else if (c[0]== CGAL::EQUAL || c[1]== CGAL::EQUAL) {
    process_aligned_centers_two_sphere_event(k,l,f, c);

  } else if (f) {
    process_intersect_event(k,l);
  } else {
    process_unintersect_event(k,l);
  }
}

void Slice::process_three_sphere_event(T::Key k, T::Key l, T::Key m, bool f) {
  if (has_degeneracy()) {
    sim_->new_event(sim_->current_time(), Three_sphere_event(this, k, l, m, f));
    process_degeneracy();
    return;
  }
  std::cout << "Three sphere kds " << k << " " << l << " " << m << " " 
	    << f << std::endl;
  Intersection_3 i3(k,l,m);
  if (!f) {
    intersections_3_.erase(i3);
  } else {
    intersections_3_[i3].first=Event_key();
  }
  CGAL_assertion(0);
}

//must maintain hafledges_

void Slice::process_rule_collapse_event(T::Key k, Rule_direction rule, T::Key o) {
  if (has_degeneracy()) {
    sim_->new_event(sim_->current_time(), Rule_event(this, k, rule, o));
    process_degeneracy();
    return;
  }
  std::cout << "Rule kds " << k << " " << rule << " " << o << std::endl;
  Rule_event_rep rep(k, rule, o);
  CGAL_precondition(rule_events_.find(rep) != rule_events_.end());
  if (rule_events_[rep].first == Event_key()) {
    rule_events_[rep].second= Event_key();
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
    uncollapse_rule(sim_->current_time(), h, base);
  } else {
    if (h->next()->next()->curve().key() == o) {
      CGAL_assertion(h->next()->next()->curve().is_arc());
      base= h->next();
    } else {
      base= h->opposite()->prev();
    }
    
    collapse_rule(sim_->current_time(), h, base);
  }
}

void Slice::process_vertex_crossing_event(Halfedge_handle h) {
  if (has_degeneracy()) {
    sim_->new_event(sim_->current_time(), Edge_event(this, h));
    process_degeneracy();
    return;
  }
  std::cout << "Rule edge "; sds_.write(h, std::cout);
  std::cout << std::endl;
  sds_.unset_event(h);

  if (h->curve().is_rule()) {
    Halfedge_handle hd= h;
    if (!h->vertex()->point().is_rule_rule()) {
      hd=h->opposite();
    }
    Halfedge_handle nh= sds_.next_edge_on_curve(hd);
    if (nh == Halfedge_handle()) {
      Halfedge_handle remove= hd->next()->opposite();
      if (remove->curve().is_outward_rule()) {
	remove= hd->opposite()->prev();
      }
      CGAL_precondition(!remove->curve().is_outward_rule());
      CGAL_precondition(remove->curve().key() != hd->curve().key());
      /*Halfedge_handle nr=*/ rotate_rule(sim_->current_time(), remove);
    }
  }

  Halfedge_handle p;
  Halfedge_handle t;
  // bool lefty;
  do {
    if (h->next()->curve().is_rule() && h->next()->curve() != h->curve()) {
      //lefty=false;
      p= h->next()->opposite();
      t= h->prev();
      break;
    } else if (h->prev()->curve().is_rule() && h->prev()->curve() != h->curve()) {
      //lefty=true;
      p=h->prev();
      t=h->next()->opposite();
      break;
    } else {
      h=h->opposite();
    } 
  } while(true);

  Halfedge_handle nt=insert_vertex(p->vertex()->point(), t);
  move_edge_target(p, nt);
  check_remove_vertex(h);

  p= sds_.cross_edge(h)->opposite();
  CGAL_assertion(p != Halfedge_handle());
  if (!p->curve().is_rule()) {
    // h=h->opposite();
    p= sds_.cross_edge(h->opposite())->opposite();
  }
  CGAL_assertion(p->curve().is_rule());
  // simplify?
  if (h->opposite()->face() == p->face() 
      || h->opposite()->face() == p->opposite()->face()) {
    h=h->opposite();
  }
  
  Qt_examiner_viewer_2 *qt= new Qt_examiner_viewer_2();
  draw_rz(qt, CGAL::to_double(sim_->current_time()) + .1);
  qt->show_everything();
  qt->show();
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
      Event_key ek1;
      if (ep.second.is_valid()) {
	// for degeneracy
	ek1= sim_->new_event(ep.second, Two_sphere_event(this,
							 i2.first(),
							 i2.second(),
							 false));
      }
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
  // There is an edge collapse if this edge and its endpoints are
  // defined by an arc and two rules (of different direction) or
  // by two arcs and one rule (h must be an arc in this case)

  // complicated a bit elimitate two SS endpoints
  std::cout << "Checking edge collapse for ";
  sds_.write(h, std::cout) << std::endl;
  CGAL_precondition(sds_.event(h)== Event_key());
  if (!h->curve().is_finite()) return;
  // an arc on a sphere with its two rules (but they might be called
  // something else
  if (h->vertex()->point().is_sphere_extremum() 
      && h->opposite()->vertex()->point().is_sphere_extremum())  return;
  
  T::Key rules[2];
  T::Key arc, oarc;
  if (h->curve().is_arc()) {
    arc= h->curve().key();
  } else {
    rules[project(h->curve().constant_coordinate())]= h->curve().key();
  }

  {
    Vertex_handle v= h->vertex();
    if (!v->point().other_curve_is_rule(h->curve())) {
      oarc= hvc.key();
    } else {
      int pc= project(v->point().other_curve_rule_coordinate(h->curve()));
      if (rules[pc] == T::Key()) {
	rules[pc]= hvc.key();
      } else {
	// two rules in the same direction
	std::cout << "Two same rules." << std::endl;
	return;
      }
    }
  }

  {
    Vertex_handle v= h->opposite()->vertex();
    if (!v->point().other_curve_is_rule(h->curve())) {
      oarc= hvc.key();
    } else {
      int pc= project(v->point().other_curve_rule_coordinate(h->curve()));
      if (rules[pc] == T::Key()) {
	rules[pc]= hvc.key();
      } else {
	// two rules in the same direction
	std::cout << "Two same rules." << std::endl;
	return;
      }
    }
  }

  T::Event_pair ep;
  if (oarc == Sds::Curve::Key()) {
    CGAL_assertion(rules[0] != T::Key() && rules[1] != T::Key());
    T::Event_pair ep= t_.sphere_intersect_rule_rule_events(arc,
							   rules[0],
							   rules[1]);
  } else {
    CGAL_assertion(rules[0] != T::Key() || rules[1] != T::Key());
    CGAL_assertion(rules[0] == T::Key() || rules[1] == T::Key());
    Coordinate_index ci;
    T::Key rule;
    if (rules[0] != T::Key()) {
      rule=rules[0];
      ci= plane_coordinate(0);
    } else {
      rule= rules[1];
      ci= plane_coordinate(1);
    }
    ep= t_.circle_cross_rule_events(arc, oarc, rule, ci);
  }
  
  Event_key ek= sim_->null_event();
  // degeneracy?
  if (ep.first.is_valid() && ep.first != ep.second) {
    if (ep.first  > sim_->current_time()) {
      ek = sim_->new_event(ep.first, Edge_event(this, h));
    } else if (ep.second > sim_->current_time()) {
      ek= sim_->new_event(ep.second, Edge_event(this,h));
    }
  }
  sds_.set_event(h, ek);
  std::cout << "Created event " << ek << std::endl;
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
      Rule_direction ri= r->curve().rule_direction();
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
    ret= p0.compare(p1, plane_coordinate(1));
    return ret;
    /*if (ret != CGAL::EQUAL) {
      return ret;
    } else {
      ret=  CGAL::compare(sim_->event<Event_base>(k0).type(),
			  sim_->event<Event_base>(k1).type());
      return ret;
      }*/
  }
}
