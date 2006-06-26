#include <CGAL/Arrangement_of_spheres_3/Slice_data_structure.h>
#include <map>
#include <set>
#include <vector>
#include <CGAL/Union_find.h>
#include <CGAL/HalfedgeDS_items_decorator.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/HalfedgeDS_const_decorator.h>

//typedef CGAL::HalfedgeDS_items_decorator<HDS> HDSID;
//typedef CGAL::HalfedgeDS_decorator<HDS> HDSD;

Slice_data_structure::Slice_data_structure() {
  initialize();
  //std::cout << hds_.size_of_bfaces() << " " << hds_.size_of_halfedges() << " " << hds_.size_of_vertices() << std::endl;
  //hds_.clear();
}

Slice_data_structure::HDS& Slice_data_structure::hds() {
  return hds_;
}

unsigned int Slice_data_structure::degree(Vertex_const_handle v) const {
  unsigned int r=0;
  Halfedge_const_handle h= v->halfedge();
  do {
    ++r;
    h= h->opposite()->prev();
  } while(h != v->halfedge());
  return r;
}


Slice_data_structure::Halfedge_handle
Slice_data_structure::split_face(Halfedge_handle o,
				 Halfedge_handle d,
				 Curve c) {
  /*Halfedge_handle h= new_halfedge(c);
  h->set_vertex(d);
  h->opposite()->set_vertex(o);
  h->set_next(d->next());
  h->opposite()->set_next(o->next());
  d->set_next(h->opposite());
  o->set_next(h);
  h->opposite()->set_face(o->face());
  o->face()->set_halfedge(h->opposite());
  Face_handle f=  Face_handle f= hds_.faces_push_back(HDS::Face());*/
  CGAL_precondition(o->face() == d->face());
  std::cout << "Spliting face ";
  write_face(o, std::cout) << std::endl;
  CGAL::HalfedgeDS_decorator<HDS> hdsd(hds_);
  Halfedge_handle h= hdsd.split_face(o,d);
  h->set_curve(c);
  h->opposite()->set_curve(c.other_side());
  std::cout << "Got ";
  write_face(h, std::cout) << " and ";
  write_face(h->opposite(), std::cout) << std::endl;
  return h;
}

Slice_data_structure::Halfedge_handle 
Slice_data_structure::remove_redundant_vertex(Vertex_handle v) {
  std::cout << "Removing redundant vertex " << v->point() 
	    << " from edges " << v->halfedge()->curve()
	    << " and " << v->halfedge()->opposite()->prev()->curve() << std::endl;
  CGAL_assertion(degree(v)==2);
  CGAL_assertion(v->halfedge()->curve().other_side()
		 == v->halfedge()->opposite()->prev()->curve());
  Halfedge_handle h= v->halfedge();
  Halfedge_handle oh= h->next();
  connect(h, oh->next());
  connect(oh->opposite()->prev(), h->opposite());
  h->face()->set_halfedge(h);
  h->opposite()->face()->set_halfedge(h->opposite());
  
  h->set_vertex(oh->vertex());
  h->vertex()->set_halfedge(h);
  hds_.edges_erase(oh);
  hds_.vertices_erase(v);
  return h;
}

Slice_data_structure::Face_handle
Slice_data_structure::remove_rule(Halfedge_handle h) {
  if (h->face() != h->opposite()->face()) {

    std::cout << "Merging faces ";
    write_face(h, std::cout) << " and ";
    write_face(h->opposite(), std::cout) << std::endl;

    Face_handle f= h->face();

    Vertex_handle s= h->vertex();
    Vertex_handle t= h->opposite()->vertex();

    CGAL::HalfedgeDS_decorator<HDS> dec(hds_);
    dec.join_face(h);

    if (degree(s) ==2 
	&& s->halfedge()->curve() == s->halfedge()->opposite()->prev()->curve()) {
      remove_redundant_vertex(s);
    }
    if (degree(t) ==2
	&& t->halfedge()->curve() == t->halfedge()->opposite()->prev()->curve()) {
      remove_redundant_vertex(t);
    }

    std::cout << "Got face ";
    write_face(f->halfedge(), std::cout) << std::endl;
    return f;
  } else {
    // check that vertex degree is 1
    std::cout << "Removing peninsula " << h->curve() << "--" 
	      << h->vertex()->point() << std::endl;
    {
      int deg=0;
      Halfedge_handle hc= h->vertex()->halfedge();
      do {
	++deg;
	h= h->opposite()->prev();
      } while (hc != h->vertex()->halfedge());
      CGAL_assertion(deg==1);
    }
    
    
    Halfedge_handle hp= h->prev();
    Halfedge_handle hn= h->next()->next();
    Face_handle f= h->face();
    Vertex_handle vg= hp->vertex();
    std::cout << "Ground is " << hp->curve() << "--" << vg->point()
	      << "--" << hn->curve() << std::endl;

    h->vertex()->point()= Point();
    CGAL_precondition(h->vertex()->halfedge()==h);
    CGAL_precondition(h->next()->opposite()==h);
    CGAL_precondition(hp->face()==f);
    CGAL_precondition(hn->face()==f);

    hds_.vertices_erase(h->vertex());
    
    CGAL_assertion(hp->vertex()==vg);
    vg->set_halfedge(hp);
    connect(hp, hn);
    
    f->set_halfedge(hp);

    hds_.edges_erase(h);
    std::cout << "Got face ";
    write_face(f->halfedge(), std::cout) << std::endl;

    {
      std::cout << "Auditing const decorator..." << std::flush;
      CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
      if (!chds.is_valid(false, 3)) {
	CGAL_assertion(0);
	std::cerr << "Not valid." << std::endl;
      }
      std::cout << "done." << std::endl;
    }
    return f;
  }
}

Slice_data_structure::Halfedge_handle 
Slice_data_structure::find_halfedge(Vertex_handle v, Face_handle f) {
  Halfedge_handle e=v->halfedge(), h=e;
  do {
    if (h->face()==f) return h;
    h= h->opposite()->prev();
  } while (h != e);
  CGAL_assertion(0);
  return Halfedge_handle();
}


Slice_data_structure::Vertex_handle 
Slice_data_structure::insert_vertex_in_edge(Halfedge_handle h, Point p) {
  Halfedge_handle i=h->opposite();

  Vertex_handle vh= h->vertex();
  Vertex_handle vi= i->vertex();
  
  Halfedge_handle hp=h->prev();
  Halfedge_handle in=i->next();

  Vertex_handle v= new_vertex(p);
  Halfedge_handle j= new_halfedge(h->curve());
  Halfedge_handle k= j->opposite();

  v->set_halfedge(i);
  vi->set_halfedge(k);
 
  j->set_face(h->face());
  k->set_face(i->face());

  k->set_vertex(vi);
  i->set_vertex(v);
  j->set_vertex(v);
  connect(j, h);
  connect(i, k);
  connect(k, in);
  connect(hp, j);
  return v;
}

void Slice_data_structure::connect(Halfedge_handle a, Halfedge_handle b) {
  a->set_next(b);
  b->set_prev(a);
}

void
Slice_data_structure::new_circle(Curve::Key k, Face_handle fi,
				 Halfedge_handle vs[4]) {
  Face_handle f= hds_.faces_push_back(HDS::Face());
  Vertex_handle vhs[4];
  vhs[0]= new_vertex(Point(Curve(k, Curve::RT_ARC),
			   Curve(k, Curve::R_RULE)));
  vhs[1]= new_vertex(Point(Curve(k, Curve::RT_ARC),
			   Curve(k, Curve::T_RULE)));
  vhs[2]= new_vertex(Point(Curve(k, Curve::LT_ARC),
			   Curve(k, Curve::L_RULE)));
  vhs[3]= new_vertex(Point(Curve(k, Curve::RB_ARC),
			   Curve(k, Curve::B_RULE)));

  vs[0]=new_halfedge(Curve(k, Curve::RT_ARC));
  vs[1]=new_halfedge(Curve(k, Curve::LT_ARC));
  vs[2]=new_halfedge(Curve(k, Curve::LB_ARC));
  vs[3]=new_halfedge(Curve(k, Curve::RB_ARC));
  
  Halfedge_handle ivs[4];
  for (unsigned int i=0; i< 4; ++i) {
    vs[i]->set_vertex(vhs[i]);
    ivs[i]= vs[i]->opposite();
    ivs[i]->set_vertex(vhs[(i+1)%4]);
    vhs[(i+1)%4]->set_halfedge(ivs[i]);
  }

  for (unsigned int i=0; i< 4; ++i) {
    connect(ivs[i], ivs[(i+1)%4]);
    ivs[i]->set_face(f);
  }

  f->set_halfedge(ivs[0]);
  CGAL_assertion(ivs[0]->next()->next()->next()->next() == ivs[0]);
  CGAL_assertion(ivs[0]->prev()->prev()->prev()->prev() == ivs[0]);

  {
    for (unsigned int i=0; i< 4; ++i) {
      connect(vs[(i+1)%4], vs[i]);
      vs[i]->set_face(fi);
    }
    fi->set_halfedge(vs[0]);
    CGAL_assertion(vs[0]->next()->next()->next()->next() == vs[0]);
    CGAL_assertion(vs[0]->prev()->prev()->prev()->prev() == vs[0]);


    std::cout << "Auditing const decorator..." << std::flush;
    CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
    if (!chds.is_valid(true, 3)) {
      std::cerr << "Not valid." << std::endl;
    }
    std::cout << "done." << std::endl;
    std::cout << "Outside: ";
    write_face(vs[0], std::cout) << std::endl;
    std::cout << "Inside: ";
    write_face(ivs[0], std::cout) << std::endl;

    {
      std::cout << "Auditing const decorator..." << std::flush;
      CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
      if (!chds.is_valid(false, 3)) {
	CGAL_assertion(0);
	std::cerr << "Not valid." << std::endl;
      }
      std::cout << "done." << std::endl;
    }
  }
}

void
Slice_data_structure::new_target(Curve::Key k, 
				 Halfedge_handle ts[4]) {

  // NOTE eventually I want to cache these
  Face_handle ft[4];
  for (unsigned int i=0; i< 4; ++i){
    ft[i]= hds_.faces_push_back(HDS::Face());
  }
  Halfedge_handle c[4];
  new_circle(k, ft[0], c);
 
  /*CGAL_assertion(tr->curve().is_top() && tr->curve().is_right()
    && !tr->curve().is_inside());
    Halfedge_handle tl= tr->opposite()->next()->opposite();
    CGAL_assertion(tl->curve().is_top() && tl->curve().is_left()
    && !tl->curve().is_inside());
    Halfedge_handle bl= tl->opposite()->next()->opposite();
    CGAL_assertion(bl->curve().is_bottom() && bl->curve().is_left()
    && !bl->curve().is_inside());
    Halfedge_handle br= bl->opposite()->next()->opposite();
    CGAL_assertion(br->curve().is_bottom() && br->curve().is_right()
    && !br->curve().is_inside());*/

  // make all point out
  ts[0]= new_halfedge(Curve::make_rule(k, 0));
  ts[1]= new_halfedge(Curve::make_rule(k, 1))->opposite();
  ts[2]= new_halfedge(Curve::make_rule(k, 2))->opposite();
  ts[3]= new_halfedge(Curve::make_rule(k, 3));
  for (unsigned int i=0; i< 4; ++i) {
    ts[i]->opposite()->set_vertex(c[i]->vertex());
    connect(c[i], ts[i]);
    connect(ts[(i+1)%4]->opposite(), c[i]);
  }

  // to make it a complete hds
  //Vertex_handle rv;
  {
    Vertex_handle rv= new_vertex(Point::make_special(k));
    rv->set_halfedge(ts[0]);
    
    for (unsigned int i=0; i< 4; ++i){
      ts[i]->set_vertex(rv);
      connect(ts[i], ts[(i+1)%4]->opposite());
    }

    CGAL::HalfedgeDS_items_decorator<HDS> dec;
    for (unsigned int i=0; i< 4; ++i){
      dec.set_face_in_face_loop(c[i], ft[i]);
      ft[i]->set_halfedge(c[i]);
    }
    write_face(c[0], std::cout) << std::endl;
    write_face(c[1], std::cout) << std::endl;
    write_face(c[2], std::cout) << std::endl;
    write_face(c[3], std::cout) << std::endl;

    std::cout << "Auditing const decorator..." << std::flush;
    CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
    if (!chds.is_valid(false, 3)) {
      CGAL_assertion(0);
      std::cerr << "Not valid." << std::endl;
    }
    std::cout << "done." << std::endl;
  }
  
  //return rv;
}

void
Slice_data_structure::insert_target(Halfedge_handle ts[4],
				    Halfedge_handle vs[4]){
  //audit(2);

  Halfedge_handle vsn[4];
  for (unsigned int i=0; i< 4; ++i){
    CGAL_assertion(vs[i]->face() == vs[(i+1)%4]->face());
    vsn[i]= vs[i]->next();
  }
  Face_handle old_face= vs[0]->face();
  
  std::cout << "Inserting into ";
  write_face(vs[0], std::cout) << std::endl;

  for (unsigned int i=0; i< 4; ++i){
    std::cout << i << " is " << vs[i]->curve() << "--" << vs[i]->vertex()->point() 
	      << "--" << vsn[i]->curve() << std::endl;
  }

  Vertex_handle tempv= ts[0]->vertex();
  
  for (unsigned int i=0; i< 4; ++i){
    ts[i]->set_vertex(vs[i]->vertex());
    connect(ts[i], vsn[i]);
    connect(vs[i], ts[i]->opposite());
   
  }

  hds_.vertices_erase(tempv);
  
  std::cout << "Got:" << std::endl;
  for (unsigned int i=0; i< 4; ++i) {
    write_face(vs[i], std::cout) << std::endl;
  }

  CGAL::HalfedgeDS_items_decorator<HDS> dec;
  for (unsigned int i=0; i< 4; ++i) {
    Face_handle f = ts[i]->face();
   
    dec.set_face_in_face_loop(vsn[i], f);
  }

  hds_.faces_erase(old_face);
  
  std::cout << "Auditing structure..." << std::flush;
  audit(1); 
  std::cout << "done." << std::endl;
}


Slice_data_structure::Face_handle 
Slice_data_structure::remove_target(Halfedge_handle ts[4]) {
  audit(1);
  // NOTE later cache this

  Halfedge_handle out;
  for (unsigned int i=0; i< 4; ++i) {
    std::cout << ts[i]->curve() << "--" << ts[i]->vertex()->point() << std::endl;
    connect(ts[i]->opposite()->prev(), ts[i]->next());
    out=ts[i]->next();
    ts[i]->vertex()->set_halfedge(ts[i]->opposite()->prev());
  
  }
  for (unsigned int i=0; i< 4; ++i) {
    if (degree(ts[i]->vertex())==2) {
      out=remove_redundant_vertex(ts[i]->vertex());
    }
  }

  Face_handle f= hds_.faces_push_back(HDS::Face());
  f->set_halfedge(out);
  CGAL::HalfedgeDS_items_decorator<HDS> dec;
  dec.set_face_in_face_loop(out, f);
  std::cout << "Got face ";
  write_face(out, std::cout) << std::endl;

  {
    Vertex_handle rv= new_vertex(Point::make_special(ts[0]->curve().key()));
    for (unsigned int i=0; i< 4; ++i){
      connect(ts[i], ts[(i+1)%4]->opposite());
      ts[i]->set_vertex(rv);
      ts[i]->face()->set_halfedge(ts[i]);
    }
    rv->set_halfedge(ts[0]);
  }
  
 
  
 

  audit(2);
  {
    // destroy target
    hds_.faces_erase(ts[0]->prev()->opposite()->face());
    hds_.vertices_erase(ts[0]->vertex());
    
    for (unsigned int i=0; i< 4; ++i) {
      hds_.vertices_erase(ts[i]->opposite()->vertex());
      hds_.faces_erase(ts[i]->face());
      hds_.edges_erase(ts[i]->prev());
      hds_.edges_erase(ts[i]);
    }
  }
  audit(1);
  return f;
}

bool Slice_data_structure::has_vertex(Face_const_handle fh, 
				      Vertex_const_handle vh) const {
  Halfedge_const_handle h= vh->halfedge();
  do {
    if (h->face() == fh) return true;
    h= h->opposite()->prev();
    CGAL_assertion(h->vertex() == vh);
  } while (h != vh->halfedge());
  return false;
}


void Slice_data_structure::audit_halfedge(Halfedge_const_handle h) const {
  CGAL_assertion(h->next() != Halfedge_const_handle());
  CGAL_assertion(h->prev() != Halfedge_const_handle());
  CGAL_assertion(h->opposite() != Halfedge_const_handle());
  CGAL_assertion(h->next()->prev()==h);
  CGAL_assertion(h->prev()->next()==h);
  CGAL_assertion(h->opposite()->opposite()==h);
  CGAL_assertion(h->face() != Face_const_handle());
  CGAL_assertion(h->vertex() != Vertex_const_handle());
  CGAL_assertion(h->opposite()->vertex() != h->vertex());
  CGAL_assertion(h->opposite()->face() != h->face());
}

void Slice_data_structure::audit_vertex(Vertex_const_handle v, bool has_special) const {
  Point pt = v->point();
  
  if (pt.is_special()) {
    CGAL_assertion(has_special);
    std::cout << "Skipping special point." << std::endl;
    return;
  }
  //  bool printed_point=false;
  // check for tangency points
  std::vector<Curve> curves;
  HDS::Halfedge_const_handle c= v->halfedge();
  do {
    curves.push_back(c->curve());
    c= c->next()->opposite();
  } while (c != v->halfedge());
    
  CGAL_assertion(!curves.empty());
  CGAL_assertion(curves.size() != 1);
  if (curves.size() ==2) {
    if (/*curves[0]== curves[1] ||*/ curves[0].other_side() == curves[1]){
      std::cerr << "Warning, collinear vertex in " << curves[0] << std::endl;
    } else {
      CGAL_assertion(!curves[0].is_finite());
      CGAL_assertion(!curves[1].is_finite());
      CGAL_assertion(curves[0].is_rule());
      CGAL_assertion(curves[1].is_rule());
      CGAL_assertion(curves[0].is_vertical()
		     != curves[1].is_vertical());
    }
  } else {
      
    //CGAL_assertion(curves.size() ==4);
    std::map<Curve::Key,int> arcs;
    std::map<Curve::Key,std::pair<int,bool> > rules;
    std::vector<Curve> ordered_arcs;
    bool fa=false, fb=false;
    for (unsigned int i=0; i< curves.size(); ++i){
      Curve c= curves[i];
      //      Curve nc= curves[(i+1)%curves.size()];
      Curve::Key ind= c.key();
      if (curves[i]== pt.first() || curves[i].other_side() == pt.first()){
	fa=true;
      }
      if (curves[i]== pt.second() || curves[i].other_side() == pt.second()){
	fb=true;
      }
      //Curve::Key nind= nc.key();
      /*if (!(c.is_rule() || nc.is_rule() || ind != nind)){
	std::cerr << "Vertex " << v->point() << " is improperly defined";
	}*/
      if (c.is_rule()){
	if (rules.find(ind)== rules.end()){
	  rules[ind]= std::make_pair(0, c.is_vertical());
	} else {
	  CGAL_assertion(rules[ind].first==0);
	  ++rules[ind].first;
	  CGAL_assertion(rules[ind].second== c.is_vertical());
	}
      } else {
	if (arcs.find(ind) == arcs.end()){
	  arcs[ind]=1;
	} else {
	  ++arcs[ind];
	}
	ordered_arcs.push_back(c);
      }
    }
    //CGAL_assertion(fa);
    //CGAL_assertion(fb);
    if (!fa || !fb) {
      std::cerr << "Vertex " << pt << " is surrounded by edges ";
      for (unsigned int i=0; i< curves.size(); ++i){
	std::cerr << curves[i] << "--";
      }
      std::cerr << std::endl;
      CGAL_assertion(0);
    }
    for (std::map<Curve::Key,int>::const_iterator it= arcs.begin();
	 it != arcs.end(); ++it){
      CGAL_assertion(it->second ==2);
    }
      
      
    // check that all arc pairs cross
    CGAL_assertion(ordered_arcs.size()%2==0);
    int half= ordered_arcs.size()/2;
    for (unsigned int i=0; i< ordered_arcs.size(); ++i){
      CGAL_assertion(ordered_arcs[i].key() 
		     == ordered_arcs[(i+half)%ordered_arcs.size()].key());
    }
  }
}



struct Handle_compare{
  template <class H>
  bool operator()(H a, H b) const {
    return &*a < &*b;
  }
};

void Slice_data_structure::audit(unsigned int num_set) const {

  std::set<HE_key> reachable;
  std::cout << "Auditing vertices..." << std::flush;
  for (HDS::Vertex_const_iterator it= hds_.vertices_begin(); 
       it != hds_.vertices_end(); ++it){

    audit_vertex(it, num_set != 1);

    HDS::Halfedge_const_handle c= it->halfedge();
    int ct=0;
    do {
      reachable.insert(c);
      ++ct;
      c= c->next()->opposite();
      if (ct > 100) {
	std::cerr << "Taking too long to walk around vertex " 
		  << it->point() << std::endl;
	int nct=0;
	HDS::Halfedge_const_handle nc= it->halfedge();
	do {
	  ++nct;
	  std::cout << nc->curve() << std::endl;
	  if (nct==10) break;
	  nc= nc->next()->opposite();
	} while (nc != it->halfedge());
	CGAL_assertion(0);
      }
    } while (c != it->halfedge());
  }
  std::cout << "done." << std::endl;

  std::cout << "Auditing halfedges..." << std::flush;
  for (HDS::Halfedge_iterator it= hds_.halfedges_begin();
       it != hds_.halfedges_end(); ++it){
    HDS::Halfedge_const_handle h= it;
    HDS::Halfedge_const_handle ho= it->opposite();

    if (it->next() == HDS::Halfedge_handle()) {
      std::cerr<< "Invalid next for ";
      write(it, std::cerr) << std::endl;
      errors_.push_back(it->curve());
    }
    if (it->prev() == HDS::Halfedge_handle()) {
      std::cerr<< "Invalid prev for ";
      write(it, std::cerr) << std::endl;
      errors_.push_back(it->curve());
    } else if (it != it->prev()->next()){
      std::cerr<< "Invalid prev/next for ";
      write(it, std::cerr) << std::endl;
      errors_.push_back(it->curve());
    }
    if (reachable.find(h) == reachable.end()){
      std::cerr << "Non-vertex reachable halfedge ";
      write(it, std::cerr) << std::endl;
      errors_.push_back(it->curve());
    }
      
    CGAL_assertion(it->curve().is_valid());
    CGAL_assertion(it->curve().is_inside() != it->opposite()->curve().is_inside());
  }

  std::cout << "done." << std::endl;
  
  std::cout << "Auditing faces..." << std::flush;
  for (HDS::Face_iterator it= hds_.faces_begin();
       it != hds_.faces_end(); ++it){
    Face_handle f=it;
    CGAL_assertion(f->halfedge() != Halfedge_handle());
  }
  std::cout << "done." << std::endl;

  std::cout << "Auditing connectivity..." << std::flush;
  CGAL::Union_find<Vertex_const_handle> uf;
  std::map<Vertex_const_handle, 
    CGAL::Union_find<Vertex_const_handle>::handle,
    Handle_compare > handles;

  for (HDS::Vertex_const_iterator it= hds_.vertices_begin(); 
       it != hds_.vertices_end(); ++it){
    handles[it]= uf.make_set(it);
  }


  for (HDS::Halfedge_iterator it= hds_.halfedges_begin();
       it != hds_.halfedges_end(); ++it){
    uf.unify_sets(handles[it->vertex()],
		  handles[it->opposite()->vertex()]);
  }
  std::cout << "done." << std::endl;

 
  CGAL_assertion(uf.number_of_sets()==num_set);

  std::cout << "Auditing const decorator..." << std::flush;
  CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
  if (!chds.is_valid(true, 3)) {
    CGAL_assertion(0);
    std::cerr << "Not valid." << std::endl;
  }
  std::cout << "done." << std::endl;
    
}

void Slice_data_structure::set_is_building(bool tf) {
  if (tf==false) {
    for (std::map<Edge, Halfedge_handle>::iterator it= unmatched_hedges_.begin(); it != unmatched_hedges_.end(); ++it){
      //std::cout << "Searching for next for ";
      //write(it->second, std::cout) << std::endl;
      it->second->set_face(inf_);
      inf_->set_halfedge(it->second);
      Halfedge_handle c=it->second->opposite()->prev()->opposite();
      Vertex_handle v= it->second->vertex();
      while (c->prev() != Halfedge_handle()){
	//write( c, std::cout) << std::endl;
	Vertex_handle vo= c->opposite()->vertex();
	CGAL_assertion(v==vo);
	c= c->prev()->opposite();
      }
      //std::cout << "Found ";
      //write(c, std::cout) << std::endl;
      c->set_prev(it->second);
      it->second->set_next(c);
    }
    unmatched_hedges_.clear();
  }
}
  
std::ostream &Slice_data_structure::write(Halfedge_const_handle h, 
					  std::ostream &out) const {
  out << h->opposite()->vertex()->point() << " -- " << h->curve()
      << " -- " << h->vertex()->point();
  return out;
}

std::ostream &Slice_data_structure::write_face(Halfedge_const_handle h,
					       std::ostream &out) const {
  Halfedge_const_handle e= h;
  do {
    out << h->curve() << "--" << h->vertex()->point() << "--" << std::flush;
    h=h->next();
  } while (h != e);
  return out;
}

void Slice_data_structure::reserve(int nv, int ne, int nf) {
  hds_.reserve(nv, ne, nf);
}

void Slice_data_structure::initialize() {
  inf_= hds_.faces_push_back(HDS::Face());
}

void Slice_data_structure::clear() {
  hds_.vertices_clear();
  hds_.edges_clear();
  hds_.faces_clear();
  inf_=Face_handle();
  errors_.clear();
    
  unmatched_hedges_.clear();
  points_.clear();
  initialize();
}

//typedef std::pair<Point,Point>  ED;

Slice_data_structure::Vertex_handle Slice_data_structure::new_vertex(Point p) {
  CGAL_precondition(p.is_valid());
  //std::cout << "Creating point " << p << std::endl;
  HDS::Vertex v;
  v.point()=p;
  points_[p]=hds_.vertices_push_back(v);
  return points_[p];
}

Slice_data_structure::Halfedge_handle Slice_data_structure::new_halfedge( Curve ff){
  Slice_data_structure::Halfedge_handle h=hds_.edges_push_back(HDS::Halfedge(), HDS::Halfedge());
  h->set_curve(ff);
  h->opposite()->set_curve(ff.other_side());
  return h;
}

Slice_data_structure::Halfedge_handle Slice_data_structure::new_halfedge(Point s, Curve ff, Point f) {
  //std::cout << "Creating edge " << s << " -- " << ff << " -- " << f << std::endl;
  CGAL_precondition(points_.find(s) != points_.end());
  if (points_.find(f) == points_.end()) {
    new_vertex(f);
  }

  Halfedge_handle h;
  Edge ep(s, ff, f);
  if (unmatched_hedges_.find(ep) != unmatched_hedges_.end()){
    h= unmatched_hedges_[ep];
    unmatched_hedges_.erase(ep);
    CGAL_assertion(h->is_border());
    //std::cout << "matched" << std::endl;
  } else {
    //std::cout << "unmatched" << std::endl;
    h= new_halfedge(ff);
    //h->set_inside(inside);
    //h->opposite()->set_inside(!inside);
    unmatched_hedges_[Edge(f, ff.other_side(), s)]= h->opposite();
  }
  points_[f]->set_halfedge(h);
  h->set_vertex(points_[f]);
  points_[s]->set_halfedge(h->opposite());
  h->opposite()->set_vertex(points_[s]);
  return h;
}


