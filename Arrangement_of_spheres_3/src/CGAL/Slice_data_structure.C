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

Slice_data_structure::Slice_data_structure(int num) {
  initialize(num);
  //std::cout << hds_.size_of_bfaces() << " " << hds_.size_of_halfedges() << " " << hds_.size_of_vertices() << std::endl;
  //hds_.clear();
}

Slice_data_structure::HDS& Slice_data_structure::hds() {
  return hds_;
}

Slice_data_structure::Halfedge_handle 
Slice_data_structure::a_halfedge(Curve::Key k) const {
  CGAL_precondition(halfedges_[k.input_index()][0]->curve().is_inside());
  return halfedges_[k.input_index()][0];
}

Slice_data_structure::Halfedge_handle 
Slice_data_structure::rule_halfedge(Curve::Key k, int i) const {
  Halfedge_handle h= halfedges_[k.input_index()][i]->opposite()->prev();
  CGAL_assertion(h->vertex()->point().is_sphere_extremum());
  CGAL_assertion(h->curve().key() == k);
  CGAL_assertion(h->curve().is_rule());
  return h;
}

Slice_data_structure::Halfedge_handle 
Slice_data_structure::next_edge_on_curve(Halfedge_handle h) const {
  if (h->curve().is_rule()) {
    if (!h->vertex()->point().is_rule_rule()) return Halfedge_handle();
    int deg = degree(h->vertex());
    Halfedge_handle r= h->next();
    do {
      if (// check if the rules go in the same direction
	  r->curve().constant_coordinate() == h->curve().constant_coordinate()
	  // if the degree is 4 we can stop if they are not the same key
	  // otherwise we have to continue in certain cases
	  && (deg ==3 || r->curve().key() == h->curve().key())) {
	return r;
      }
      r=r->opposite()->next();
    } while (r != h->opposite());
    return Halfedge_handle();
  } else {
    Halfedge_handle r= h->next();
    do {
      if (h->curve().is_rule() == r->curve().is_rule()
	  && h->curve().key() == r->curve().key()) {
	return r;
      }
      r=r->opposite()->next();
    } while (r != h->opposite());
    return Halfedge_handle();
  }
}

Slice_data_structure::Halfedge_const_handle 
Slice_data_structure::next_edge_on_curve(Halfedge_const_handle h) const {
  CGAL_assertion(0);
  return h;
}

Slice_data_structure::Halfedge_handle 
Slice_data_structure::cross_edge(Halfedge_handle h) const {
  CGAL_precondition(degree(h->vertex()) == 3 || degree(h->vertex()) ==4 );
  Halfedge_handle n= h->next();
  if (h->curve().is_rule() && n->curve().is_rule() 
      && h->curve().constant_coordinate() == n->curve().constant_coordinate()
      || !h->curve().is_rule() && n->curve().key() ==h->curve().key() 
      && n->opposite()->curve().is_rule() ==h->curve().is_rule()) {
    n= h->opposite()->prev()->opposite();
  } 
  CGAL_assertion(n->curve() !=h->curve() && n->opposite()->curve() != h->curve());
  return n;
}

bool Slice_data_structure::is_redundant(Vertex_const_handle v) const {
  if (degree(v) != 2) return false;
  if (v->halfedge()->curve().is_rule()) {
    if (v->halfedge()->curve().constant_coordinate() 
	!= v->halfedge()->next()->curve().constant_coordinate()) return false;
    else return true;
  } else {
    return v->halfedge()->curve() == v->halfedge()->next()->curve();
  }
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
Slice_data_structure::remove_redundant_vertex(Halfedge_handle h) {
  std::cout << "Removing redundant vertex " << h->vertex()->point() 
	    << " from edges " << h->curve()
	    << " and " << h->opposite()->prev()->curve() << std::endl;
  Vertex_handle v= h->vertex();
  CGAL_assertion(degree(v)==2);
  if (h->curve().is_rule()) {
    CGAL_assertion(h->curve().is_vertical()
		   == h->opposite()->prev()->curve().is_vertical());
  } else {
    CGAL_assertion(h->curve().other_side()
		   == h->opposite()->prev()->curve());
  }
  Halfedge_handle nh= h->next();
  connect(h->prev(), nh);
  connect(nh->opposite(), h->opposite()->next());
  h->face()->set_halfedge(nh);
  h->opposite()->face()->set_halfedge(nh->opposite());
  
  nh->opposite()->set_vertex(h->opposite()->vertex());
  nh->opposite()->vertex()->set_halfedge(nh->opposite());
  hds_.edges_erase(h);
  hds_.vertices_erase(v);
  return nh;
}

std::pair<Slice_data_structure::Halfedge_handle,
	  Slice_data_structure::Halfedge_handle>
Slice_data_structure::remove_rule(Halfedge_handle h) {
  if (h->face() != h->opposite()->face()) {

    std::cout << "Merging faces ";
    write_face(h, std::cout) << " and ";
    write_face(h->opposite(), std::cout) << std::endl;

    Face_handle f= h->face();

    Vertex_handle s= h->vertex();
    Vertex_handle t= h->opposite()->vertex();
    Halfedge_handle sh= h->opposite()->prev();
    Halfedge_handle th= h->prev();

    CGAL::HalfedgeDS_decorator<HDS> dec(hds_);
    h->set_curve(Curve());
    dec.join_face(h);

    
    /*if (degree(s) ==2) {
      CGAL_assertion(s->halfedge()->curve()
		     == s->halfedge()->opposite()->prev()->opposite()->curve());
      // fails for degeneracies, I'll fix that elsewhere (in the degeneracy handler)
      remove_redundant_vertex(s->halfedge());
    } else {
      std::cout << "Not removing vertex " <<s->point() << std::endl;
    }
    if (degree(t) ==2) {
      CGAL_assertion(t->halfedge()->curve() 
		     == t->halfedge()->opposite()->prev()->opposite()->curve());
      remove_redundant_vertex(t->halfedge());
    } else {
      std::cout << "Not removing vertex " << t->point() << std::endl;
      }*/

    std::cout << "Got face ";
    write_face(f->halfedge(), std::cout) << std::endl;
    std::cout << "And edge ";
    write(sh, std::cout);
    std::cout << " and ";
    write(th, std::cout);
    std::cout << std::endl;
    return std::make_pair(sh, th);
  } else {
    
    // check that vertex degree is 1
    std::cout << "Removing peninsula " << h->curve() << "--" 
	      << h->vertex()->point() << std::endl;
    CGAL_assertion(0);
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
    return std::make_pair(Halfedge_handle(), Halfedge_handle());
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
  for (unsigned int i=0; i< 4; ++i) {
    vhs[i]= new_vertex(Point::make_extremum(k, i));
  }

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
    if (!chds.is_valid(false, 3)) {
      std::cerr << "Not valid." << std::endl;
    }
    std::cout << "done." << std::endl;
    std::cout << "Outside: ";
    write_face(vs[0], std::cout) << std::endl;
    std::cout << "Inside: ";
    write_face(ivs[0], std::cout) << std::endl;

    if (0) {
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

  if (targets_.empty()) {
    
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
  } else {
    // recycle
    Vertex_handle rv= targets_.back();
    targets_.pop_back();
    Halfedge_handle h= rv->halfedge();
    do {
      ts[h->curve().rule_index()]= h;
      h= h->next()->opposite();
    } while (h != rv->halfedge());
  
    relabel_target(ts, k);
  }
  
  for (unsigned int i=0; i< 4; ++i) {
    CGAL_assertion(ts[i]->event() == Simulator::Event_key());
    CGAL_assertion(ts[i]->prev()->event() == Simulator::Event_key());
  }
  //return rv;
}

void
Slice_data_structure::insert_target(Curve::Key k,
				    Halfedge_handle vs[4]){
  //audit(2);

  Halfedge_handle ts[4];
  new_target(k, ts);

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
  
  for (unsigned int i=0; i< 4; ++i) {
    halfedges_[k.input_index()][i]= vs[i]->next()->next()->opposite();
  }
  
  std::cout << "Auditing structure..." << std::flush;
  audit(); 
  std::cout << "done." << std::endl;
}

bool Slice_data_structure::is_in_slice(Vertex_const_handle v) const{
  return !v->point().is_special();
}

bool Slice_data_structure::is_in_slice(Halfedge_const_handle h) const{
  return !h->curve().key().is_target();
}

bool Slice_data_structure::is_in_slice(Face_const_handle h) const{
  return is_in_slice(h->halfedge());
}


void
Slice_data_structure::relabel_target(Halfedge_handle ts[], Curve::Key k) {
  // 8 edges, 4 vertices
  CGAL_precondition(k.is_target() || ts[0]->curve().key().is_target());
  for (unsigned int i=0; i< 4; ++i) {
    //std::cout << "Relabeling " << ts[i]->curve() << " to " << k << std::endl;
    ts[i]->curve().set_key(k);
    //std::cout << "Got " << ts[i]->curve() << " with " << k << std::endl;
    CGAL_assertion(ts[i]->curve().key() == k );

    //std::cout << "Relabeling " << ts[i]->opposite()->curve() << " to " << k << std::endl;
    ts[i]->opposite()->curve().set_key(k);
    ts[i]->opposite()->vertex()->point().set_key(k);
    //std::cout << "Relabeling " << ts[i]->prev()->curve() << " to " << k << std::endl;
    ts[i]->prev()->curve().set_key(k);
    //std::cout << "Relabeling " << ts[i]->prev()->opposite()->curve() 
    //	      << " to " << k << std::endl;
    ts[i]->prev()->opposite()->curve().set_key(k);

    CGAL_assertion(ts[i]->curve().key() == k );
  }
  CGAL_assertion(ts[0]->curve().key() == k );
  CGAL_assertion(ts[1]->curve().key() == k );
  CGAL_assertion(ts[2]->curve().key() == k );
  CGAL_assertion(ts[3]->curve().key() == k );
}

Slice_data_structure::Face_handle 
Slice_data_structure::remove_target(Halfedge_handle ts[4],
				    Halfedge_handle verts[4]) {
  audit();
  Curve::Key k= ts[0]->curve().key();
  // NOTE later cache this

  Halfedge_handle out;
  for (unsigned int i=0; i< 4; ++i) {
    std::cout << ts[i]->curve() << "--" << ts[i]->vertex()->point() << std::endl;
    connect(ts[i]->opposite()->prev(), ts[i]->next());
    out=ts[i]->next();
    ts[i]->vertex()->set_halfedge(ts[i]->opposite()->prev());
    verts[i]= ts[i]->opposite()->prev();
  }
  /*for (unsigned int i=0; i< 4; ++i) {
    if (degree(ts[i]->vertex())==2) {
      out=remove_redundant_vertex(ts[i]->vertex());
    }
    }*/

  Face_handle f= hds_.faces_push_back(HDS::Face());
  f->set_halfedge(out);
  CGAL::HalfedgeDS_items_decorator<HDS> dec;
  dec.set_face_in_face_loop(out, f);
  std::cout << "Got face ";
  write_face(out, std::cout) << std::endl;

  {
    Vertex_handle rv= new_vertex(Point::make_special(Curve::Key::target_key()));
    for (unsigned int i=0; i< 4; ++i){
      connect(ts[i], ts[(i+1)%4]->opposite());
      ts[i]->set_vertex(rv);
      ts[i]->face()->set_halfedge(ts[i]);
    }
    rv->set_halfedge(ts[0]);
    targets_.push_back(rv);
    relabel_target(ts, Curve::Key::target_key());
    CGAL_assertion(ts[0]->curve().key().is_target());
    CGAL_assertion(ts[1]->curve().key().is_target());
    CGAL_assertion(ts[2]->curve().key().is_target());
    CGAL_assertion(ts[3]->curve().key().is_target());
    CGAL_assertion(rv->halfedge()->curve().key().is_target());
  }
  
 
  
 

 
  /*{
    // destroy target
    hds_.faces_erase(ts[0]->prev()->opposite()->face());
    hds_.vertices_erase(ts[0]->vertex());
    
    for (unsigned int i=0; i< 4; ++i) {
      hds_.vertices_erase(ts[i]->opposite()->vertex());
      hds_.faces_erase(ts[i]->face());
      hds_.edges_erase(ts[i]->prev());
      hds_.edges_erase(ts[i]);
    }
    }*/
  halfedges_[k.input_index()][0]= Halfedge_handle();
  halfedges_[k.input_index()][1]= Halfedge_handle();
  halfedges_[k.input_index()][2]= Halfedge_handle();
  halfedges_[k.input_index()][3]= Halfedge_handle();
  

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


/*void Slice_data_structure::audit_halfedge(Halfedge_const_handle h) const {

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
  }*/

void Slice_data_structure::audit_vertex(Vertex_const_handle v) const {
  Point pt = v->point();
  
  if (!is_in_slice(v)) {
    /*
      Special vertices are, for example, ones in targets that are not
      in use.
    */
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
    
  if (v->point().is_sphere_extremum()) {
    CGAL_assertion(degree(v) ==3 || degree(v)==4);
  }

  CGAL_assertion(!curves.empty());
  CGAL_assertion(curves.size() != 1);
  if (curves.size() ==2) {
    if (/*curves[0]== curves[1] ||*/ curves[0].other_side() == curves[1]){
      /* There should never be a rule with an unused vertex in it */
      std::cerr << "Warning, collinear vertex in " << curves[0] << std::endl;
      std::cerr << v->point() << std::endl;
      std::cerr << curves[0] << std::endl;
      std::cerr << curves[1] << std::endl;
      CGAL_assertion(0);
    } else {
      /*
	The only vertices of degree 2 should be the four corners. 
      */
      CGAL_assertion(!curves[0].is_finite());
      CGAL_assertion(!curves[1].is_finite());
      CGAL_assertion(curves[0].is_rule());
      CGAL_assertion(curves[1].is_rule());
      CGAL_assertion(curves[0].is_vertical()
		     != curves[1].is_vertical());
    }
  } else {
    std::map<Curve::Key,int> arcs;
    std::map<Curve::Key,std::pair<int,bool> > rules;
    std::vector<Curve> ordered_arcs;

    /*
      There can only be one of each rule entering a vertex.
      Each arc that enters a vertex leaves the vertex.
    */
    for (unsigned int i=0; i< curves.size(); ++i){
      Curve c= curves[i];
      Curve::Key ind= c.key();

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

void Slice_data_structure::audit() const {
  unsigned int num_set= 1+ targets_.size();

  std::cout << "Auditing const decorator..." << std::flush;
  CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
  if (!chds.is_valid(true, 3)) {
    CGAL_assertion(0);
    std::cerr << "Not valid." << std::endl;
  }
  std::cout << "done." << std::endl;


  std::set<HE_key> reachable;
  std::cout << "Auditing vertices..." << std::flush;
  for (HDS::Vertex_const_iterator it= hds_.vertices_begin(); 
       it != hds_.vertices_end(); ++it){

    audit_vertex(it);

    HDS::Halfedge_const_handle c= it->halfedge();
    int ct=0;
    do {
      reachable.insert(c);
      ++ct;
      c= c->next()->opposite();
      // produce better error messages
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

    /*
      All halfedges are connected to some valid vertex.
    */
    if (reachable.find(h) == reachable.end()){
      std::cerr << "Non-vertex reachable halfedge ";
      write(it, std::cerr) << std::endl;
      errors_.push_back(it->curve());
    }
      
    CGAL_assertion(it->curve().key().is_target() || it->curve().is_valid());
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

  /* The number of connected componenets is 1+ # of saved targets*/

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

  /*
    Everything in a saved target is marked as a target and 
    nothing edges is.
  */

  std::cout << "Auditing targets..." << std::flush;
  for (HDS::Halfedge_iterator it= hds_.halfedges_begin();
       it != hds_.halfedges_end(); ++it){
    if (it->curve().key().is_target()) {
      CGAL_assertion(it->vertex()->point().is_special());
      CGAL_assertion(it->vertex()->point().is_special());
      CGAL_assertion(it->opposite()->curve().key().is_target());
    } else {
      CGAL_assertion(!it->vertex()->point().is_special());
      CGAL_assertion(!it->vertex()->point().is_special());
      CGAL_assertion(!it->opposite()->curve().key().is_target());
    }
  }

  for (unsigned int i=0; i< targets_.size(); ++i){
    CGAL_assertion(targets_[i]->point().is_special());
  }
  std::cout << "done." << std::endl;
  
  std::cout << "Auditing halfedge map..." << std::flush;
  //CGAL_assertion(halfedges_.size() == t_.number_of_spheres());
  for (unsigned int i=0; i< halfedges_.size(); ++i){
    for (unsigned int j=0; j<4; ++j) {
      if (halfedges_[i][0] != Halfedge_handle()) {
	CGAL_assertion(halfedges_[i][j]->curve().key() == Curve::Key(i));
	CGAL_assertion(halfedges_[i][j]->curve().is_inside());
	CGAL_assertion(halfedges_[i][j]->opposite()->prev()->curve().is_rule());
	//CGAL_assertion(halfedges_[i][j]->opposite()->prev()->curve().key()==T::Key(i));
      } else {
	CGAL_assertion(halfedges_[i][j] == Halfedge_handle());
      }
    }
  }

  for (Halfedge_const_iterator vit= halfedges_begin();
       vit != halfedges_end(); ++vit) {
    if (vit->curve().key() != Curve::Key::target_key()) {
       if (vit->curve().is_arc()) {
	 CGAL_assertion(halfedges_[vit->curve().key().input_index()][0]
			!= Halfedge_handle());
       }
     }
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

    for (Halfedge_iterator it= halfedges_begin(); 
       it != halfedges_end(); ++it){
    if (it->curve().is_arc() && it->curve().key().is_input()
	&& it->curve().is_inside()) {
      Halfedge_handle fit= next_edge_on_curve(it);
      if (fit->curve() != it->curve()) {
	int ai= (it->curve().arc_index()+1)%4;
	halfedges_[it->curve().key().input_index()][ai]=it;
      }
    }
  }

  } else {
    hds_.vertices_clear();
    hds_.edges_clear();
    hds_.faces_clear();
    inf_= hds_.faces_push_back(HDS::Face());
    errors_.clear();
    
    unmatched_hedges_.clear();
    points_.clear();
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

void Slice_data_structure::initialize(int num) {
  inf_= hds_.faces_push_back(HDS::Face());
  Face_handle in= hds_.faces_push_back(HDS::Face());
  Curve tc(Point::Key(Point::Key::BL), Curve::T_RULE);
  Curve rc(Point::Key(Point::Key::BL), Curve::R_RULE);
  Curve bc(Point::Key(Point::Key::TR), Curve::B_RULE);
  Curve lc(Point::Key(Point::Key::TR), Curve::L_RULE);
  Halfedge_handle lh= new_halfedge(lc);
  Halfedge_handle bh= new_halfedge(bc);
  Halfedge_handle th= new_halfedge(tc);
  Halfedge_handle rh= new_halfedge(rc);
  Point blp(rc, tc);
  Point tlp(lc, tc);
  Point trp(lc, bc);
  Point brp(rc, bc);
  Vertex_handle blv= new_vertex(blp);
  Vertex_handle tlv= new_vertex(tlp);
  Vertex_handle trv= new_vertex(trp);
  Vertex_handle brv= new_vertex(brp);
  
  bh->set_vertex(brv);
  brv->set_halfedge(bh);
  rh->set_vertex(brv);
  
  bh->opposite()->set_vertex(trv);
  trv->set_halfedge(bh->opposite());
  lh->set_vertex(trv);

  lh->opposite()->set_vertex(tlv);
  tlv->set_halfedge(lh->opposite());
  th->opposite()->set_vertex(tlv);

  th->set_vertex(blv);
  blv->set_halfedge(th);
  rh->opposite()->set_vertex(blv);

  connect(lh->opposite(), th);
  connect(th, rh);
  connect(rh, bh->opposite());
  connect(bh->opposite(), lh->opposite());
  connect(bh, rh->opposite());
  connect(rh->opposite(), th->opposite());
  connect(th->opposite(), lh);
  connect(lh, bh);

  CGAL::HalfedgeDS_items_decorator<HDS> dec;
  dec.set_face_in_face_loop(lh, in);
  in->set_halfedge(lh);
  CGAL_assertion(lh->opposite()->face() == Face_handle());
  dec.set_face_in_face_loop(lh->opposite(), inf_);
  inf_->set_halfedge(lh->opposite());
  //audit();
  halfedges_.resize(num);
}

void Slice_data_structure::clear() {
  hds_.vertices_clear();
  hds_.edges_clear();
  hds_.faces_clear();
  inf_=Face_handle();
  errors_.clear();
    
  unmatched_hedges_.clear();
  points_.clear();
  initialize(halfedges_.size());
}

//typedef std::pair<Point,Point>  ED;

Slice_data_structure::Vertex_handle Slice_data_structure::new_vertex(Point p) {
  //CGAL_precondition(p.is_valid() || p.is_special());
  //std::cout << "Creating point " << p << std::endl;
  Vertex_handle v= hds_.vertices_push_back(HDS::Vertex());
  v->point()=p;
  return v;
}

Slice_data_structure::Vertex_handle Slice_data_structure::new_vertex_cached(Point p) {
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
    new_vertex_cached(f);
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


void Slice_data_structure::exchange_vertices(Halfedge_handle h,
					     Halfedge_handle p){

  std::cout << "Exchanging vertices with ";
  write(h, std::cout);
  std::cout << " and ";
  write(p, std::cout);
  std::cout << " as input" << std::endl;
  CGAL_precondition(!p->vertex()->point().is_sphere_extremum());
  if (p->vertex() == h->vertex()) {
    Halfedge_handle n= p->next();
    Halfedge_handle t= h->prev();
    Halfedge_handle f= t->opposite()->prev();
    Halfedge_handle g= h->opposite()->next();
    CGAL_precondition(h->curve() == n->curve());
    CGAL_precondition(p->curve().is_rule());
    
    connect(f, h->opposite());
    connect(h->opposite(), t->opposite());
    connect(t, p->opposite());
    connect(p, h);
    connect(h, n);
    connect(n->opposite(), g);
    
    h->set_face(p->face());
    h->opposite()->set_face(f->face());
    
    Vertex_handle va= p->vertex();
    Vertex_handle vb= t->vertex();
    //f->set_vertex(va);
    CGAL_assertion(f->vertex() == vb);
    t->set_vertex(va);
    h->set_vertex(vb);
    va->set_halfedge(t);
    vb->set_halfedge(h);
    va->point()= Point(p->curve(), t->curve());

    h->opposite()->set_vertex(va);
    n->opposite()->set_vertex(vb);

    h->set_curve(t->curve());
    h->opposite()->set_curve(t->opposite()->curve());
  } else {
    CGAL_precondition(p->vertex() == h->opposite()->vertex());
    CGAL_precondition(p==h->prev());
    //Halfedge_handle p= h->prev();
    Halfedge_handle n= h->next();
    Halfedge_handle t= p->opposite()->prev();
    Halfedge_handle f= h->opposite()->prev();
    Halfedge_handle g= n->opposite()->next();
    CGAL_precondition(h->curve() == t->curve());
    CGAL_precondition(p->curve().is_rule());
    
    connect(f, t->opposite());
    connect(t, h);
    connect(h, p->opposite());
    connect(p, n);
    connect(n->opposite(), h->opposite());
    connect(h->opposite(), g);
    
    h->opposite()->set_face(g->face());
    h->set_face(t->face());

    Vertex_handle va= t->vertex();
    Vertex_handle vb= h->vertex();
    va->point()= Point(p->curve(), n->curve());
    //f->set_vertex(va);
    CGAL_assertion(f->vertex() == vb);
    t->set_vertex(vb);
    h->set_vertex(va);
    vb->set_halfedge(t);
    va->set_halfedge(h);
    h->opposite()->set_vertex(vb);
    n->opposite()->set_vertex(va);

    h->set_curve(n->curve());
    h->opposite()->set_curve(n->opposite()->curve());
  }
}



Slice_data_structure::Face_handle
Slice_data_structure::intersect(Halfedge_handle ha, Halfedge_handle hb) {
  CGAL_precondition(ha->face() == hb->face());
  Halfedge_handle ha_prev= ha->prev();
  Halfedge_handle ha_op_next= ha->opposite()->next();
  Halfedge_handle ha_next= ha->next();
  Halfedge_handle ha_op_prev= ha->opposite()->prev();

  Halfedge_handle hb_prev= hb->prev();
  Halfedge_handle hb_op_next= hb->opposite()->next();
  Halfedge_handle hb_next= hb->next();
  Halfedge_handle hb_op_prev= hb->opposite()->prev();
  
  Vertex_handle va= new_vertex(Point(ha->curve(), hb->curve()));
  Vertex_handle vb= new_vertex(Point(hb->curve(), ha->curve()));

  Halfedge_handle hap= new_halfedge(ha->curve());
  Halfedge_handle hbp= new_halfedge(hb->curve());
  Halfedge_handle fha= new_halfedge(ha->curve());
  Halfedge_handle fhb= new_halfedge(hb->curve());

  Face_handle fn= hds_.faces_push_back(HDS::Face());
  Face_handle fa= hds_.faces_push_back(HDS::Face());
  Face_handle f= ha->face();
  
  Vertex_handle ha_v= ha->vertex();
  Vertex_handle hb_v= hb->vertex();
  Vertex_handle ha_op_v= ha->opposite()->vertex();
  Vertex_handle hb_op_v= hb->opposite()->vertex();

  hap->set_vertex(ha_v);
  hap->opposite()->set_vertex(vb);
  hb->set_vertex(vb);
  hb->opposite()->set_vertex(hb_op_v);

  fha->set_vertex(vb);
  fha->opposite()->set_vertex(va);

  fhb->opposite()->set_vertex(vb);
  fhb->set_vertex(va);

  hbp->opposite()->set_vertex(va);
  ha->set_vertex(va);
  ha->opposite()->set_vertex(ha_op_v);
  hbp->set_vertex(hb_v);

  va->set_halfedge(fhb);
  vb->set_halfedge(fha);
  ha_v->set_halfedge(hap);
  hb_op_v->set_halfedge(hb->opposite());
  ha_op_v->set_halfedge(ha->opposite());
  hb_v->set_halfedge(hbp);
  

  connect(ha, hbp);
  connect(hbp, hb_next);
  connect(hb_op_prev, hbp->opposite());
  connect(hbp->opposite(), fha);
  connect(fhb, ha->opposite());
  connect(fha->opposite(), fhb->opposite());
  connect(fhb->opposite(), fha->opposite());
  connect(hap->opposite(), fhb);
  connect(hb, hap);
  connect(fha, hb->opposite());
  connect(hap, ha_next);
  connect(ha_op_prev, hap->opposite());

 

  fhb->set_face(ha->opposite()->face());
  fhb->opposite()->set_face(fn);
  fha->set_face(hb->opposite()->face());
  fha->opposite()->set_face(fn);
  hap->set_face(f);
  hap->opposite()->set_face(ha->opposite()->face());
  hbp->set_face(fa);
  hbp->opposite()->set_face(fha->face());
  ha->set_face(fa);

  fa->set_halfedge(ha);
  CGAL::HalfedgeDS_items_decorator<HDS> dec;
  dec.set_face_in_face_loop(ha, fa);

  fhb->face()->set_halfedge(fhb);
  fha->face()->set_halfedge(fha);
  fn->set_halfedge(fha->opposite());
  f->set_halfedge(hb);

  write_face(fn->halfedge(), std::cout);
  std::cout << std::endl;
  write_face(fha->face()->halfedge(), std::cout);
  std::cout << std::endl;
  write_face(fhb->face()->halfedge(), std::cout);
  std::cout << std::endl;
  write_face(ha->face()->halfedge(), std::cout);
  std::cout << std::endl;
  write_face(hb->face()->halfedge(), std::cout);
  std::cout << std::endl;

  return fn;
}


std::pair<Slice_data_structure::Halfedge_handle,Slice_data_structure::Halfedge_handle>
Slice_data_structure::unintersect(Face_handle fn) {
  Halfedge_handle fha= fn->halfedge()->opposite();
  Halfedge_handle fhb= fha->opposite()->next()->opposite();
  CGAL_assertion(fhb->opposite()->next()->opposite() == fha);

  Vertex_handle va= fhb->vertex();
  Vertex_handle vb= fha->vertex();
  Halfedge_handle hap= fhb->prev()->opposite();
  Halfedge_handle hb= fha->next()->opposite();
  Halfedge_handle hbp= fha->prev()->opposite();
  Halfedge_handle ha = fhb->next()->opposite();
  Vertex_handle ha_v= hap->vertex();
  //Vertex_handle hb_op_v= hb->opposite()->vertex();
  Vertex_handle hb_v = hbp->vertex();
  Face_handle f= hb->face();
  Face_handle fa= ha->face();
  
  Halfedge_handle ha_next= hap->next();
  Halfedge_handle ha_op_prev= hap->opposite()->prev();
  Halfedge_handle hb_next= hbp->next();
  Halfedge_handle hb_op_prev = hbp->opposite()->prev();

  CGAL::HalfedgeDS_items_decorator<HDS> dec;
  dec.set_face_in_face_loop(ha, f);

  // delete things
  hds_.faces_erase(fn);
  hds_.faces_erase(fa);
  hds_.edges_erase(fhb);
  hds_.edges_erase(fha);
  hds_.edges_erase(hap);
  hds_.edges_erase(hbp);
  hds_.vertices_erase(va);
  hds_.vertices_erase(vb);

  connect(ha, ha_next);
  connect(ha_op_prev, ha->opposite());
  connect(hb, hb_next);
  connect(hb_op_prev, hb->opposite());
  ha->set_face(f);
  f->set_halfedge(ha);
  ha->opposite()->face()->set_halfedge(ha->opposite());
  hb->opposite()->face()->set_halfedge(hb->opposite());
  
  ha->set_vertex(ha_v);
  ha_v->set_halfedge(ha);
  hb->set_vertex(hb_v);
  hb_v->set_halfedge(hb);
  return std::make_pair(ha, hb);
}


void Slice_data_structure::move_rule(Halfedge_handle r, 
				     Halfedge_handle ne, 
				     Halfedge_handle t) {
  // remove r
  // stitch across base
  // insert new base and vertex
  // attach r
  CGAL_precondition(r->vertex()->point().is_sphere_extremum());
  CGAL_precondition(r->curve().is_rule());

  Vertex_handle vr= r->vertex();
  Curve cr= r->curve();
  //near= t->next()->next();
  
  Halfedge_handle nv= insert_vertex_in_edge(ne, vr->point())->halfedge();
  
  CGAL::HalfedgeDS_decorator<HDS> hdsd(hds_);
  hdsd.join_face(r);
  r= hdsd.split_face(t, nv);
  r->set_curve(cr);
  r->opposite()->set_curve(cr.other_side());
  remove_redundant_vertex(vr->halfedge());
}
