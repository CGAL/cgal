#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>

#include <map>
#include <set>
#include <vector>
#include <CGAL/Union_find.h>
#include <CGAL/HalfedgeDS_items_decorator.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/HalfedgeDS_const_decorator.h>
#include <CGAL/Arrangement_of_spheres_3/Rule_direction.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE


//typedef CGAL::HalfedgeDS_items_decorator<HDS> HDSID;
//typedef CGAL::HalfedgeDS_decorator<HDS> HDSD;

CGAL_AOS3_TEMPLATE
Combinatorial_cross_section CGAL_AOS3_TARG::Combinatorial_cross_section() {
  //std::cout << hds_.size_of_bfaces() << " " << hds_.size_of_halfedges() << " " << hds_.size_of_vertices() << std::endl;
  //hds_.clear();
  clear();
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::HDS&
Combinatorial_cross_section CGAL_AOS3_TARG::hds() {
  return hds_;
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle 
Combinatorial_cross_section CGAL_AOS3_TARG::a_halfedge(Curve::Key k) const {
  CGAL_precondition(halfedges_[k.input_index()]== Halfedge_handle()
		    || halfedges_[k.input_index()]->curve().is_inside() 
		    && halfedges_[k.input_index()]->curve().key() == k);
  return halfedges_[k.input_index()];
}


CGAL_AOS3_TEMPLATE
void 
Combinatorial_cross_section CGAL_AOS3_TARG::set_halfedge(Halfedge_handle h){
  halfedges_[h->curve().key().input_index()] = h;
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle 
Combinatorial_cross_section CGAL_AOS3_TARG::next_edge_on_rule(Halfedge_handle h) const {
  CGAL_precondition(h->curve().is_rule());
  if (!h->vertex()->point().is_rule_rule()) return Halfedge_handle();
  //int deg = degree(h->vertex());
  Halfedge_handle r= h->next();
  do {
    if (r->curve().constant_coordinate() == h->curve().constant_coordinate()) {
      return r;
    }
    r=r->opposite()->next();
  } while (r != h->opposite());
  return Halfedge_handle();
}



CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle 
Combinatorial_cross_section CGAL_AOS3_TARG::next_edge_on_curve(Halfedge_handle h) const {
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
      if (h->curve().is_arc() == r->curve().is_arc()
	  && h->curve().key() == r->curve().key()) {
	CGAL_assertion(h->curve().is_inside()
		       == r->curve().is_inside());
	return r;
      }
      r=r->opposite()->next();
      CGAL_assertion(r!= h->opposite());
    } while (true);
    CGAL_assertion(0);
    return Halfedge_handle();
  }
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle 
Combinatorial_cross_section CGAL_AOS3_TARG::cross_edge(Halfedge_handle h) const {
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

CGAL_AOS3_TEMPLATE
bool Combinatorial_cross_section CGAL_AOS3_TARG::is_redundant(Vertex_const_handle v) const {
  if (degree(v) != 2) return false;
  if (v->halfedge()->curve().is_rule()) {
    if (v->halfedge()->curve().constant_coordinate() 
	!= v->halfedge()->next()->curve().constant_coordinate()) return false;
    else return true;
  } else {
    return v->halfedge()->curve() == v->halfedge()->next()->curve();
  }
}

CGAL_AOS3_TEMPLATE
unsigned int Combinatorial_cross_section CGAL_AOS3_TARG::degree(Vertex_const_handle v) const {
  unsigned int r=0;
  Halfedge_const_handle h= v->halfedge();
  do {
    ++r;
    h= h->opposite()->prev();
  } while(h != v->halfedge());
  return r;
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle
Combinatorial_cross_section CGAL_AOS3_TARG::split_face(Curve c,
				 Halfedge_handle o,
				 Halfedge_handle d){
  
  CGAL_precondition(o->face() == d->face());
  std::cout << "Spliting face ";
  write(o->face(), std::cout) << std::endl;
  CGAL::HalfedgeDS_decorator<HDS> hdsd(hds_);
  Halfedge_handle h= hdsd.split_face(o,d);
  set_curve(h, c);
  std::cout << "Got ";
  write(h->face(), std::cout) << " and ";
  write(h->opposite()->face(), std::cout) << std::endl;
  return h;
}



CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle 
Combinatorial_cross_section CGAL_AOS3_TARG::remove_vertex(Halfedge_handle h) {
  std::cout << "Removing vertex " << h->vertex()->point() 
	    << " from edges " << h->curve()
	    << " and " << h->opposite()->prev()->curve() << std::endl;
  Vertex_handle v= h->vertex();
  CGAL_assertion(degree(v)==2);
  //if (labels_ok) {
    if (h->curve().is_rule()) {
      CGAL_assertion(h->curve().is_vertical()
		     == h->opposite()->prev()->curve().is_vertical());
    } else {
      CGAL_assertion(h->curve().other_side()
		     == h->opposite()->prev()->curve());
    }
    //}
  Halfedge_handle nh= h->next();
  connect(h->prev(), nh);
  connect(nh->opposite(), h->opposite()->next());
  h->face()->set_halfedge(nh);
  h->opposite()->face()->set_halfedge(nh->opposite());
  
  nh->opposite()->set_vertex(h->opposite()->vertex());
  nh->opposite()->vertex()->set_halfedge(nh->opposite());
  write(nh, std::cout) << std::endl;
  write(nh->opposite(), std::cout) << std::endl;
  if (nh->opposite()->vertex()->point().is_sphere_extremum() 
	     && nh->curve().is_arc() && !nh->curve().is_inside()) {
    halfedges_[nh->curve().key().input_index()] =(nh->opposite());
  }
  h->set_vertex(Vertex_handle());
  delete_edge(h);
  hds_.vertices_erase(v);
  return nh;
}

/*std::pair<Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle,
  Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle>*/
CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Face_handle
Combinatorial_cross_section CGAL_AOS3_TARG::merge_faces(Halfedge_handle h) {
  std::cout << "Merging faces ";
  write(h->face(), std::cout) << " and ";
  write(h->opposite()->face(), std::cout) << std::endl;

  //CGAL_precondition(degree(h->vertex()) >=3);
  //CGAL_precondition(degree(h->opposite()->vertex()) >=3);

  Face_handle f= h->face();

  Vertex_handle s= h->vertex();
  Vertex_handle t= h->opposite()->vertex();
  Halfedge_handle sh= h->opposite()->prev();
  Halfedge_handle th= h->prev();
  
  CGAL::HalfedgeDS_decorator<HDS> dec(hds_);
  h->set_curve(Curve());
  dec.join_face(h);
  
  std::cout << "Got face ";
  write(f->halfedge()->face(), std::cout) << std::endl;
  std::cout << "And edge ";
  write(sh, std::cout);
  std::cout << " and ";
  write(th, std::cout);
  std::cout << std::endl;
  //turn std::make_pair(sh, th);
  return f;
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Face_handle
Combinatorial_cross_section CGAL_AOS3_TARG::merge_faces(Vertex_handle v) {
  {
    std::cout << "Merging faces ";
    Halfedge_handle h=v->halfedge();
    do {
      write(h->face(), std::cout) << std::endl;
      h= h->opposite()->prev();
    } while (h != v->halfedge()); 
  }

  Halfedge_handle h= v->halfedge();
  Face_handle f= h->face();

  
  CGAL::HalfedgeDS_decorator<HDS> dec(hds_);

  dec.erase_center_vertex(h);
  
  std::cout << "Got face ";
  write(f, std::cout) << std::endl;
  //turn std::make_pair(sh, th);
  return f;
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle 
Combinatorial_cross_section CGAL_AOS3_TARG::find_halfedge(Vertex_handle v, Face_handle f) {
  Halfedge_handle e=v->halfedge(), h=e;
  do {
    if (h->face()==f) return h;
    h= h->opposite()->prev();
  } while (h != e);
  CGAL_assertion(0);
  return Halfedge_handle();
}



CGAL_AOS3_TEMPLATE
void 
Combinatorial_cross_section CGAL_AOS3_TARG::relabel_rule(Halfedge_handle h, Curve c) {
  std::cout << "Relabeling from " << h->curve() << " to " << c << std::endl;
  //clean_edge(h);
  CGAL_precondition(h->curve().is_inside()==c.is_inside());
  CGAL_precondition(h->curve().is_vertical()==c.is_vertical());
  h->set_curve(c);
  h->opposite()->set_curve(c.other_side());
  Halfedge_handle n;
  if (h->next()->curve().is_rule() 
      && h->curve().is_vertical() == h->next()->curve().is_vertical()) {
    n= h->next();
  } else if (h->next()->curve().is_rule() 
	     && h->next()->opposite()->next()->curve().is_rule()
	     && h->curve().is_vertical() 
	     == h->next()->opposite()->next()->curve().is_vertical()) {
    n= h->next()->opposite()->next();
  }
  if (n != Halfedge_handle()) {
    CGAL_assertion(n->opposite()->vertex() == h->vertex());
    if (n->curve().is_same_part(c)) {
      relabel_rule(n, c);
    } else {
      std::cout << "Next is other dir " << std::endl; 
    }
  } else {
    std::cout << "No next " << std::endl;
  }
}





CGAL_AOS3_TEMPLATE
void
Combinatorial_cross_section CGAL_AOS3_TARG::exchange_sphere_extremums(Curve::Key k, Curve::Key l) {
  // k has the smaller center so it must be outside
  // switch them
 

  std::swap(halfedges_[k.input_index()],
	    halfedges_[l.input_index()]);
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle 
Combinatorial_cross_section CGAL_AOS3_TARG::insert_vertex(Point p, Halfedge_handle h) {
 

  Vertex_handle v= new_vertex(p);
#ifndef NDEBUG
  Face_handle f= h->face();
#endif
  Halfedge_handle r= insert_vertex_in_edge_unsafe(h,v);
  CGAL_postcondition(f == r->face());
  return r;
  //return v;
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle 
Combinatorial_cross_section CGAL_AOS3_TARG::insert_vertex_in_edge_unsafe(Halfedge_handle h,
						   Vertex_handle v) {
  Halfedge_handle i=h->opposite();

  Vertex_handle vh= h->vertex();
  Vertex_handle vi= i->vertex();
  
  Halfedge_handle hp=h->prev();
  Halfedge_handle in=i->next();

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

  if (k->curve().is_arc() && k->vertex()->point().is_sphere_extremum()
      && k->curve().is_inside()) {
    halfedges_[k->curve().key().input_index()]= k;
    std::cout << "Updated sphere halfedge to ";
    write(k, std::cout) << std::endl;
  } else {
    std::cout << "Not updating sphere halfedge ";
    write(h, std::cout) << std::endl;
    write(i, std::cout) << std::endl;
    write(j, std::cout) << std::endl;
    write(k, std::cout) << std::endl;
  }

  return j;
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME std::pair<CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle,
			     CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle>
Combinatorial_cross_section CGAL_AOS3_TARG::pinch_bl(Halfedge_handle a, Halfedge_handle b, Point p) {
  // make more efficient later
  Vertex_handle v= new_vertex(p);
  Halfedge_handle ha= insert_vertex_in_edge_unsafe(a,v);
  Halfedge_handle hb= insert_vertex_in_edge_unsafe(b,v);
  CGAL_assertion(ha->vertex() == v);
  CGAL_assertion(hb->vertex() == v);
  //CGAL_assertion(ha->face() == hb->face());
  //Vertex_handle vh= insert_vertex_in_edge(h, v->vertex()->point());
  // now merge the two vertices.
  merge_vertices_bl(ha, hb);

  std::cout << "Auditing const decorator..." << std::flush;
  CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
  if (!chds.is_valid(true, 3)) {
    CGAL_assertion(0);
    std::cerr << "Not valid." << std::endl;
  }
  std::cout << "done." << std::endl;

  return std::make_pair(ha, hb);
}


CGAL_AOS3_TEMPLATE
void Combinatorial_cross_section CGAL_AOS3_TARG::merge_vertices_bl(Halfedge_handle a,
					     Halfedge_handle b) {
  CGAL_precondition(a->face() == b->face());
  CGAL_precondition(a->vertex()->point() == b->vertex()->point());
  Vertex_handle v= a->vertex();
  Vertex_handle ov=b->vertex();
  b->set_vertex(v);
  Halfedge_handle an= a->next();
  Halfedge_handle bn= b->next();
  bn->opposite()->set_vertex(v);
  connect(b, an);
  connect(a, bn);
  Face_handle fh= hds_.faces_push_back(HDS::Face());
  fh->set_halfedge(b);
  CGAL::HalfedgeDS_items_decorator<HDS> dec;
  dec.set_face_in_face_loop(b, fh);
  if (v != ov) hds_.vertices_erase(ov);
  a->face()->set_halfedge(a);
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Vertex_handle
Combinatorial_cross_section CGAL_AOS3_TARG::unpinch_bl(Halfedge_handle a, Halfedge_handle b) {  
  CGAL_precondition(a->vertex()== b->vertex());
  CGAL_precondition(degree(a->vertex()) == 4);
  a->vertex()->set_halfedge(a);
  Halfedge_handle an= a->opposite()->prev()->opposite();
  Halfedge_handle bn= b->opposite()->prev()->opposite();
  Vertex_handle nv= new_vertex(a->vertex()->point());
  connect(a,an);
  connect(b, bn);
  nv->set_halfedge(b);
  b->set_vertex(nv);
  bn->opposite()->set_vertex(nv);
  Face_handle df= b->face();
  CGAL::HalfedgeDS_items_decorator<HDS> dec;
  dec.set_face_in_face_loop(a, a->face());
  hds_.faces_erase(df);

 std::cout << "Auditing const decorator..." << std::flush;
  CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
  if (!chds.is_valid(true, 3)) {
    CGAL_assertion(0);
    std::cerr << "Not valid." << std::endl;
  }
  std::cout << "done." << std::endl;

  return nv;
}

CGAL_AOS3_TEMPLATE
void Combinatorial_cross_section CGAL_AOS3_TARG::connect(Halfedge_handle a, Halfedge_handle b) {
  std::cout << "connecting ";
  write(a, std::cout) << " to ";
  write(b, std::cout) << std::endl;
  a->set_next(b);
  b->set_prev(a);
}

CGAL_AOS3_TEMPLATE
void
Combinatorial_cross_section CGAL_AOS3_TARG::new_circle(Curve::Key k,
						       Halfedge_handle vs[4]) {
  //audit(true);
  ++num_components_;
  /*  std::cout << "Auditing const decorator..." << std::flush;
  CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
  if (!chds.is_valid(true, 3)) {
    std::cerr << "Not valid." << std::endl;
    CGAL_assertion(0);
    }*/


  Face_handle f= hds_.faces_push_back(HDS::Face());
  Vertex_handle vhs[4];
  for (unsigned int i=0; i< 4; ++i) {
    vhs[i]= new_vertex(Point::make_extremum(k, Rule_direction(i)));
  }

  vs[0]=new_halfedge(Curve(k, Curve::RT_ARC));
  vs[1]=new_halfedge(Curve(k, Curve::LT_ARC));
  vs[2]=new_halfedge(Curve(k, Curve::LB_ARC));
  vs[3]=new_halfedge(Curve(k, Curve::RB_ARC));
  
  if (halfedges_.size() <=static_cast<unsigned int>(k.input_index())) halfedges_.resize(k.input_index()+1);
  halfedges_[k.input_index()]=vs[3]->opposite();
 
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
      vs[i]->set_face(inf_);
    }
    //fi->set_halfedge(vs[0]);
    CGAL_assertion(vs[0]->next()->next()->next()->next() == vs[0]);
    CGAL_assertion(vs[0]->prev()->prev()->prev()->prev() == vs[0]);


    /*std::cout << "Auditing const decorator..." << std::flush;
    CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
    if (!chds.is_valid(true, 3)) {
      std::cerr << "Not valid." << std::endl;
      CGAL_assertion(0);
    }
    std::cout << "done." << std::endl;*/
    std::cout << "Outside: ";
    write(vs[0]->face(), std::cout) << std::endl;
    std::cout << "Inside: ";
    write(ivs[0]->face(), std::cout) << std::endl;
  }
}

#if 0
CGAL_AOS3_TEMPLATE
void
Combinatorial_cross_section CGAL_AOS3_TARG::new_target(Curve::Key k, 
						       Halfedge_handle ts[4]) {
  //audit();
  if (targets_.empty()) {
    
    // NOTE eventually I want to cache these
    Face_handle ft[4];
    for (unsigned int i=0; i< 4; ++i){
      ft[i]= hds_.faces_push_back(CGAL_AOS3_TYPENAME HDS::Face());
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
    ts[0]= new_halfedge(Curve::make_rule(k, Rule_direction::right()));
    ts[1]= new_halfedge(Curve::make_rule(k, Rule_direction::top()))->opposite();
    ts[2]= new_halfedge(Curve::make_rule(k, Rule_direction::left()))->opposite();
    ts[3]= new_halfedge(Curve::make_rule(k, Rule_direction::bottom()));
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
      write(c[0]->face(), std::cout) << std::endl;
      write(c[1]->face(), std::cout) << std::endl;
      write(c[2]->face(), std::cout) << std::endl;
      write(c[3]->face(), std::cout) << std::endl;
      //write(c[3]->face(), std::cout) << std::endl;
      
      std::cout << "Auditing const decorator..." << std::flush;
      CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
      if (!chds.is_valid(true, 3)) {
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
      ts[h->curve().rule_direction().index()]= h;
      h= h->next()->opposite();
    } while (h != rv->halfedge());
  
    relabel_target(ts, k);
  }
  
  for (unsigned int i=0; i< 4; ++i) {
    CGAL_assertion(ts[i]->event() == Event_key());
    CGAL_assertion(ts[i]->prev()->event() == Event_key());
  }
  //return rv;
  // can't audit here since there may be collinear vertices
}


CGAL_AOS3_TEMPLATE
void
Combinatorial_cross_section CGAL_AOS3_TARG::insert_target(Curve::Key k,
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
  write(vs[0]->face(), std::cout) << std::endl;

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
    write(vs[i]->face(), std::cout) << std::endl;
  }

  CGAL::HalfedgeDS_items_decorator<HDS> dec;
  for (unsigned int i=0; i< 4; ++i) {
    Face_handle f = ts[i]->face();
   
    dec.set_face_in_face_loop(vsn[i], f);
  }

  hds_.faces_erase(old_face);
  
  for (unsigned int i=0; i< 4; ++i) {
    Halfedge_handle h= vs[i]->next()->next()->opposite();
    CGAL_assertion(h->vertex()->point().is_sphere_extremum());
    CGAL_assertion(h->vertex()->point().sphere_extremum_index() 
		   == Rule_direction(i));
    set_extremum_halfedge( h);
  }
  
  std::cout << "Auditing structure..." << std::flush;
  audit(); 
  std::cout << "done." << std::endl;
}
#endif

CGAL_AOS3_TEMPLATE
bool Combinatorial_cross_section CGAL_AOS3_TARG::is_in_slice(Vertex_const_handle v) const{
  return !v->point().is_special();
}

CGAL_AOS3_TEMPLATE
bool Combinatorial_cross_section CGAL_AOS3_TARG::is_in_slice(Halfedge_const_handle h) const{
  return h->face() != inf_; // really could check if I am the outside edge of the box
}


CGAL_AOS3_TEMPLATE
bool Combinatorial_cross_section CGAL_AOS3_TARG::is_in_slice(Face_const_handle h) const{
  return h != inf_ && is_in_slice(h->halfedge());
}


#if 0
CGAL_AOS3_TEMPLATE
void
Combinatorial_cross_section CGAL_AOS3_TARG::relabel_target(Halfedge_handle ts[], Curve::Key k) {
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


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Face_handle 
Combinatorial_cross_section CGAL_AOS3_TARG::remove_target(Halfedge_handle ts[4],
				    Halfedge_handle verts[4]) {
  audit();
  Curve::Key k= ts[0]->curve().key();
  

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
  write(out->face(), std::cout) << std::endl;

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
  halfedges_[k.input_index()]= Halfedge_handle();


  return f;
}
#endif


CGAL_AOS3_TEMPLATE
bool Combinatorial_cross_section CGAL_AOS3_TARG::has_vertex(Face_const_handle fh, 
				      Vertex_const_handle vh) const {
  Halfedge_const_handle h= vh->halfedge();
  do {
    if (h->face() == fh) return true;
    //if (h->opposite()->face() == fh) return true;
    h= h->opposite()->prev();
    CGAL_assertion(h->vertex() == vh);
  } while (h != vh->halfedge());
  return false;
}


/*void Combinatorial_cross_section CGAL_AOS3_TARG::audit_halfedge(Halfedge_const_handle h) const {

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
CGAL_AOS3_TEMPLATE
void Combinatorial_cross_section CGAL_AOS3_TARG::audit_vertex(Vertex_const_handle v, bool extra) const {
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
  CGAL_AOS3_TYPENAME HDS::Halfedge_const_handle c= v->halfedge();
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
      CGAL_assertion(extra);
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

CGAL_AOS3_TEMPLATE
void Combinatorial_cross_section CGAL_AOS3_TARG::audit(bool extra_vertices) const {
  unsigned int num_set=num_components_; //; = 1+ targets_.size();

  write(std::cout);
  std::cout << "Auditing const decorator..." << std::flush;
  CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
  if (!chds.is_valid(false, num_components_==1? 3: 2)) {
    CGAL_assertion(0);
    std::cerr << "Not valid." << std::endl;
  }
  std::cout << "done." << std::endl;


  std::set<HE_key> reachable;
  std::cout << "Auditing vertices..." << std::flush;
  for (CGAL_AOS3_TYPENAME HDS::Vertex_const_iterator it= hds_.vertices_begin(); 
       it != hds_.vertices_end(); ++it){

    audit_vertex(it, extra_vertices);

    CGAL_AOS3_TYPENAME HDS::Halfedge_const_handle c= it->halfedge();
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
	CGAL_AOS3_TYPENAME HDS::Halfedge_const_handle nc= it->halfedge();
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
  for (CGAL_AOS3_TYPENAME HDS::Halfedge_iterator it= hds_.halfedges_begin();
       it != hds_.halfedges_end(); ++it){
    CGAL_AOS3_TYPENAME HDS::Halfedge_const_handle h= it;
    CGAL_AOS3_TYPENAME HDS::Halfedge_const_handle ho= it->opposite();

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

    if (it->curve().is_arc()) {
      CGAL_assertion(next_edge_on_curve(it) != Halfedge_handle());
      CGAL_assertion(next_edge_on_curve(it->opposite()) != Halfedge_handle());
    }
      
    CGAL_assertion(it->curve().key().is_target() || it->curve().is_valid());
    CGAL_assertion(it->curve().is_inside() != it->opposite()->curve().is_inside());
  }

  std::cout << "done." << std::endl;
  
  std::cout << "Auditing faces..." << std::flush;
  for (CGAL_AOS3_TYPENAME HDS::Face_iterator it= hds_.faces_begin();
       it != hds_.faces_end(); ++it){
    Face_handle f=it;
    CGAL_assertion( f->halfedge() != Halfedge_handle());
  }
  std::cout << "done." << std::endl;

  /* The number of connected componenets is 1+ # of saved targets*/

  std::cout << "Auditing connectivity..." << std::flush;
  CGAL::Union_find<Vertex_const_handle> uf;
  CGAL_AOS3_TYPENAME  std::map<Vertex_const_handle, 
    CGAL_AOS3_TYPENAME CGAL::Union_find<Vertex_const_handle>::handle,
    Handle_compare > handles;

  for (CGAL_AOS3_TYPENAME HDS::Vertex_const_iterator it= hds_.vertices_begin(); 
       it != hds_.vertices_end(); ++it){
    handles[it]= uf.make_set(it);
  }


  for (CGAL_AOS3_TYPENAME HDS::Halfedge_iterator it= hds_.halfedges_begin();
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

  /*std::cout << "Auditing targets..." << std::flush;
  for (CGAL_AOS3_TYPENAME HDS::Halfedge_iterator it= hds_.halfedges_begin();
       it != hds_.halfedges_end(); ++it){
    if (it->curve().key().is_target()) {
      CGAL_assertion(it->vertex()->point().is_special());
      CGAL_assertion(it->vertex()->point().is_special());
      CGAL_assertion(it->opposite()->curve().key().is_target());
      CGAL_assertion(it->event() == Event_key());
    } else {
      CGAL_assertion(!it->vertex()->point().is_special());
      CGAL_assertion(!it->vertex()->point().is_special());
      CGAL_assertion(!it->opposite()->curve().key().is_target());
    }
    }*/

  /*for (unsigned int i=0; i< targets_.size(); ++i){
    CGAL_assertion(targets_[i]->point().is_special());
    }*/
  std::cout << "done." << std::endl;
  
  std::cout << "Auditing halfedge map..." << std::flush;
  //CGAL_assertion(halfedges_.size() == t_.number_of_spheres());
  for (unsigned int i=0; i< halfedges_.size(); ++i){
    if (halfedges_[i] != Halfedge_handle()) {
      CGAL_assertion(halfedges_[i]->curve().key() == Curve::Key(i));
      CGAL_assertion(halfedges_[i]->curve().is_inside());
      CGAL_assertion(halfedges_[i]->vertex()->point().is_sphere_extremum());
      
      if (!halfedges_[i]->opposite()->prev()->curve().is_rule()) {
	std::cerr << "Error in sphere halfedge for " << i 
		  << " have ";
	write(halfedges_[i], std::cerr) << " not rule is ";
	write(halfedges_[i]->opposite()->prev(), std::cerr)
	  << std::endl;
	CGAL_assertion(halfedges_[i]->opposite()->prev()->curve().is_rule());
      }
      //CGAL_assertion(halfedges_[i][j]->opposite()->prev()->curve().key()==T::Key(i));
    } 
  }

  for (Halfedge_const_iterator vit= halfedges_begin();
       vit != halfedges_end(); ++vit) {
    if (vit->curve().key() != Curve::Key::target_key()) {
       if (vit->curve().is_arc()) {
	 Curve::Key k= vit->curve().key();
	 CGAL_assertion(halfedges_[k.input_index()]
			!= Halfedge_handle());
       }
     }
   }

  std::cout << "done." << std::endl;
}

CGAL_AOS3_TEMPLATE
std::ostream &Combinatorial_cross_section CGAL_AOS3_TARG::write(Vertex_const_handle h, 
					  std::ostream &out) const {
  out << h->point();
  return out;
}


CGAL_AOS3_TEMPLATE
std::ostream &Combinatorial_cross_section CGAL_AOS3_TARG::write(Halfedge_const_handle h, 
					  std::ostream &out) const {
  if (h == Halfedge_const_handle()) out << "NULL";
  else {
    out << h->opposite()->vertex()->point() << " -- " << h->curve()
	<< " -- " << h->vertex()->point();
  }
  return out;
}


CGAL_AOS3_TEMPLATE
std::ostream &Combinatorial_cross_section CGAL_AOS3_TARG::write(Face_const_handle f,
								std::ostream &out) const {
  if (f->halfedge() == Halfedge_const_handle()) out << "empty";
  else {
    Halfedge_const_handle e= f->halfedge(), h= f->halfedge();
    do {
      out << h->curve() << "--" << h->vertex()->point() << "--" << std::flush;
      h=h->next();
    } while (h != e);
  }
 return out;
 
}

CGAL_AOS3_TEMPLATE
std::ostream &Combinatorial_cross_section CGAL_AOS3_TARG::write(std::ostream &out) const {
  out << "inf=" <<  reinterpret_cast<const void*>(&*inf_) << std::endl;
  out << hds_.size_of_vertices() << " " << hds_.size_of_halfedges() 
      << " " << hds_.size_of_faces() << std::endl;
  for (CGAL_AOS3_TYPENAME Halfedge_const_iterator it= halfedges_begin(); it != halfedges_end();
       ++it) {
    write(it, out) << ": " << reinterpret_cast<const void*>(&*it->face()) << std::endl;
  }

  for (CGAL_AOS3_TYPENAME Face_const_iterator it= faces_begin(); it != faces_end();
       ++it) {
    out << reinterpret_cast<const void*>(&*it) << ": ";
    write(it, out) << std::endl;
  }

  for (CGAL_AOS3_TYPENAME Vertex_const_iterator it= vertices_begin(); it != vertices_end();
       ++it) {
    out << it->point() << ": ";
    write(it->halfedge(), out) << std::endl;
  }
  return out;
}


CGAL_AOS3_TEMPLATE
void Combinatorial_cross_section CGAL_AOS3_TARG::reserve(int nv, int ne, int nf) {
  hds_.reserve(nv, ne, nf);
}

CGAL_AOS3_TEMPLATE
bool Combinatorial_cross_section CGAL_AOS3_TARG::has_boundary() const {
  return hds_.size_of_faces() >1;
}

CGAL_AOS3_TEMPLATE
void Combinatorial_cross_section CGAL_AOS3_TARG::set_has_boundary(bool tf) {
  if (!tf) clear();
  else {
    CGAL_precondition(hds_.size_of_vertices() ==0);
    CGAL_precondition(hds_.size_of_halfedges() ==0);
    CGAL_precondition(hds_.size_of_faces() ==1);
    std::cout << "Initializing boundary\n";
    //inf_= hds_.faces_push_back(CGAL_AOS3_TYPENAME HDS::Face());
    Face_handle in= hds_.faces_push_back(CGAL_AOS3_TYPENAME HDS::Face());
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
    num_components_=1;
    audit();
  }
 }

CGAL_AOS3_TEMPLATE
void Combinatorial_cross_section CGAL_AOS3_TARG::clear() {
  hds_.vertices_clear();
  hds_.edges_clear();
  hds_.faces_clear();
  inf_=hds_.faces_push_back(CGAL_AOS3_TYPENAME HDS::Face());
  //void_=hds_.faces_push_back(CGAL_AOS3_TYPENAME HDS::Face());
  ///hds_.faces_push_back(HDS::Face());
  errors_.clear();
  num_components_=0;
  //unmatched_hedges_.clear();
  //points_.clear();
  //initialize(halfedges_.size());
}

//typedef std::pair<Point,Point>  ED;
CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Vertex_handle 
Combinatorial_cross_section CGAL_AOS3_TARG::new_vertex(Point p) {
  //CGAL_precondition(p.is_valid() || p.is_special());
  //std::cout << "Creating point " << p << std::endl;
  Vertex_handle v= hds_.vertices_push_back(CGAL_AOS3_TYPENAME  HDS::Vertex());
  v->point()=p;
  return v;
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle 
Combinatorial_cross_section CGAL_AOS3_TARG::new_halfedge( Curve ff){
  CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle h=hds_.edges_push_back(CGAL_AOS3_TYPENAME HDS::Halfedge(),
											  CGAL_AOS3_TYPENAME HDS::Halfedge());
  h->set_curve(ff);
  h->opposite()->set_curve(ff.other_side());
  return h;
}



/*void Combinatorial_cross_section CGAL_AOS3_TARG::exchange_vertices(Halfedge_handle h,
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
  }*/

CGAL_AOS3_TEMPLATE
void Combinatorial_cross_section CGAL_AOS3_TARG::set_curve(Halfedge_handle h, Curve c) {
  h->set_curve(c);
  h->opposite()->set_curve(c.other_side());
}




CGAL_AOS3_TEMPLATE
void Combinatorial_cross_section CGAL_AOS3_TARG::delete_edge(Halfedge_handle h) {
#ifndef NDEBUG
  if (h->curve().key().is_input()) {
    CGAL_assertion(halfedges_[h->curve().key().input_index()]!=h);
    CGAL_assertion(halfedges_[h->curve().key().input_index()]!=h->opposite());
  }
  h->curve() = Curve();
  h->opposite()->curve() = Curve();
  /*CGAL_assertion(h->next()== Halfedge_handle() || h->next()->prev() != h);
  CGAL_assertion(h->prev()== Halfedge_handle() || h->prev()->next() != h);
  CGAL_assertion(h->face()== Face_handle() || h->face()->halfedge() != h);
  CGAL_assertion(h->vertex()== Vertex_handle() || h->vertex()->halfedge() != h);

  CGAL_assertion(h->opposite()->next()== Halfedge_handle() || h->opposite()->next()->prev() != h->opposite());
  CGAL_assertion(h->opposite()->prev()== Halfedge_handle() || h->opposite()->prev()->next() != h->opposite());
  CGAL_assertion(h->opposite()->face()== Face_handle() || h->opposite()->face()->halfedge() != h->opposite());
  CGAL_assertion(h->opposite()->vertex()== Vertex_handle() || h->opposite()->vertex()->halfedge() != h->opposite());*/
#endif
  hds_.edges_erase(h);
}

CGAL_AOS3_TEMPLATE
void
Combinatorial_cross_section CGAL_AOS3_TARG::move_edge_target(Halfedge_handle edge, 
				       Halfedge_handle tv) {
  Halfedge_handle ot= edge->opposite()->prev();
  tv->face()->set_halfedge(tv);
  tv->next()->face()->set_halfedge(tv->next());
  connect(edge, tv->next());
  connect(tv, edge->opposite());
  edge->set_vertex(tv->vertex());
  if (edge->face() == ot->face()){
    Halfedge_handle c= edge->next();
    do {
      CGAL_precondition(c->face() != edge->face());
      c->set_face(edge->face());
      c= c->next();
    } while (c != ot) ;
  } else {
    Halfedge_handle c= edge->opposite()->prev();
    do {
      CGAL_precondition(c->face() != edge->face());
      c->set_face(edge->opposite()->face());
      c= c->prev();
    } while (c != ot) ;
  }
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME std::pair<CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle,
	  CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle>
Combinatorial_cross_section CGAL_AOS3_TARG::intersect(Halfedge_handle ha, Halfedge_handle hb) {
  Curve ca= ha->curve();
  Curve cb= hb->curve();
  std::pair<Halfedge_handle, Halfedge_handle> p0= pinch_bl(ha, hb,
							   Point(ha->curve(),
								 hb->curve()));
  std::pair<Halfedge_handle, Halfedge_handle> p1= pinch_bl(p0.second->next(),
							   p0.second,
							   Point(p0.second->curve(),
								 p0.first->curve()));
  set_curve(p1.first, p1.second->curve().other_side());
  set_curve(p1.first->next(),  p0.first->curve().other_side());
  CGAL_precondition(p1.first->curve().key() != p1.first->next()->curve().key());
  write(p0.first, std::cout) << std::endl;
  write(p0.second, std::cout) << std::endl;
  write(p1.first, std::cout) << std::endl;
  write(p1.second, std::cout) << std::endl;
  return std::make_pair(p1.first->next(), p1.first);
}




CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME std::pair<CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle,
	  CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle>
Combinatorial_cross_section CGAL_AOS3_TARG::unintersect(Face_handle fn) {
  Halfedge_handle ha= fn->halfedge();
  Halfedge_handle hb= ha->next();
  CGAL_precondition(hb->next() == ha);
  Vertex_handle hav= ha->vertex();
  Vertex_handle hbv= hb->vertex();
  Halfedge_handle up0a= ha->opposite()->next()->opposite();
  Halfedge_handle up1b= hb->opposite()->next()->opposite();
  Vertex_handle vn0= unpinch_bl(ha, up1b);
  Vertex_handle vn1= unpinch_bl(up0a, hb);
  std::swap(ha->opposite()->curve(), hb->curve());
  std::swap(ha->curve(), hb->opposite()->curve());

  Halfedge_handle han= remove_vertex(up0a);
  Halfedge_handle hbn= remove_vertex(up1b);
  Halfedge_handle ra= remove_vertex(han);
  Halfedge_handle rb= remove_vertex(hbn);
  return std::make_pair(ra,rb);
}



CGAL_AOS3_END_INTERNAL_NAMESPACE
