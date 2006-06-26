#include <CGAL/Arrangement_of_spheres_3/Slice.h>


Slice::Face_handle Slice::insert_sphere(const T::Sphere_point_3 &cp) {
  T::Key k=t_.new_sphere(cp.sphere());
  return insert_sphere(t_.sphere_events(k).first, k);
}



Slice::Halfedge_handle 
Slice::insert_sphere_find_rule_vertex(const T::Sphere_point_3 &ep, 
				      Face_handle f,
				      Sds::Curve rule) {
  Vertex_handle v;
  try {
      Halfedge_handle h= shoot_rule(ep, f, rule);
      v=sds_.insert_vertex_in_edge(h, Sds::Point(rule, h->curve()));
  } catch (On_vertex_exception e) {
    v= e.vertex_handle();
  }
  return sds_.find_halfedge(v,f);
}


Slice::T::Key Slice::roll_back_rule(const T::Sphere_point_3 &ep,
				    Halfedge_handle cur,
				    Halfedge_handle last) {
  std::cout << "Rolling back, rolling back..." << cur->curve() << std::endl;
  // check if this is the end, if not recurse
  Halfedge_handle o=cur;
  T::Key k;
  do {
    o=o->opposite()->prev();
    CGAL_assertion(o->vertex()==cur->vertex());
    if (o== cur) break;
    if (o->opposite()->curve() == cur->curve() && o->opposite() != last) {
      std::cout << "Found one to recurse on " << std::endl;
      Halfedge_handle tmp= o->next()->opposite();
      k=roll_back_rule(ep, o->opposite(), cur);
      o= tmp;
      break;
    }
  } while (true);
  
  if ( cur->vertex()->point().type() == Sds::Point::SR) {
    CGAL_assertion(!k.is_valid());
    bool degen= (cur->next()->curve() != cur->opposite()->prev()->curve()
		 && !cur->next()->curve().is_inside());
    // we hit another sphere;
    /*CGAL_assertion(cur->vertex()->point().sphere(0).key() 
      != cur->vertex()->point().rule(0).key());*/
    if (degen) {
      k=cur->next()->curve().key();
    }
  }
  if (k.is_valid()) {
    cur->curve().flip_rule(k);
    cur->opposite()->curve().flip_rule(k);
    return k;
  } else {
    if (last == Halfedge_handle()) {
      // back to the beginning
      CGAL_assertion(cur->opposite()->vertex()->point().rule(0).key()
		     == cur->opposite()->vertex()->point().sphere(0).key());
      // redundant check
      std::cout << "Done rolling back." << std::endl;
    } else {
      // we are not the first edge
      // if not add new rule
      // remove old bit
      Halfedge_handle oe;
      Face_handle f;
      Sds::Curve ec;
      if (cur->prev()->curve().key() == cur->curve().key()) {
	oe= cur->prev();
	f= cur->face();
	ec= cur->opposite()->next()->opposite()->curve();
      } else {
	oe= cur->opposite();
	f= cur->opposite()->face();
	ec= cur->prev()->curve();
      }
     
      Halfedge_handle hv= insert_sphere_find_rule_vertex(ep, f, 
							 ec);
      std::cout << "New edge supported by " << ec 
		<< " from " << oe->vertex()->point() << " to " 
		<< hv->vertex()->point() << std::endl;
      sds_.split_face(oe, hv, ec);
      Face_handle nf= sds_.remove_rule(cur);
    
    }
    return T::Key();
  }
  CGAL_assertion(0);
}

Slice::Face_handle Slice::erase_sphere(const T::Sphere_point_3 &ep,
				       T::Key k) {
  for (Sds::Halfedge_iterator fit = sds_.halfedges_begin(); fit != sds_.halfedges_end();
       ++fit) {
    if (fit->curve().is_arc() && fit->curve().is_inside() && fit->curve().key()==k) {
      return erase_sphere(ep, fit->face());
    }
  }
  CGAL_assertion(0);
  return Face_handle();
}
								   

Slice::Face_handle Slice::erase_sphere(const T::Sphere_point_3 &ep,
				       Face_handle f) {
  T::Key k= f->halfedge()->curve().key();
  // check correctness of f
  Halfedge_handle rules[4];
  Halfedge_handle h=f->halfedge();
  int deg=0;
  do {
    ++deg;
    CGAL_assertion(h->curve().key()==k);
    CGAL_assertion(h->curve().is_arc());
    CGAL_assertion(h->curve().is_inside());
    Halfedge_handle rule= h->opposite()->next();
    int index= rule->curve().rule_index();
    CGAL_assertion(rules[index]
		   == Halfedge_handle());
    rules[index] =rule;
    h= h->next();
  } while (h != f->halfedge());
  CGAL_assertion(deg==4);
  // roll in each until I have a target in a face
  
  T::Key keys[4];
  Halfedge_handle vs[4];
  for (unsigned int i=0; i< 4; ++i){
    keys[i]= roll_back_rule(ep, rules[i], Halfedge_handle());
    if (keys[i].is_valid()) {
      vs[i]= rules[i]->opposite()->prev();
      CGAL_assertion(vs[i]->vertex() == rules[i]->vertex());
    }
  }
  
 
  sds_.remove_target(rules);
  
  // I guess I don't actually have to extend the rules any more
  // it doesn't  make anything invalid and makes the structure simpler. 
  // Unless there is only halfedge between the two circles
  bool has_split=false;
  for (unsigned int i=0; i< 4; ++i){
    if (vs[i] != Halfedge_handle() && vs[i]->curve().is_arc()) {
      // I need to handle this
      Halfedge_handle hv= insert_sphere_find_rule_vertex(ep, vs[i]->face(), 
							 Sds::Curve::make_rule(keys[i], i));
      sds_.split_face(vs[i], hv, Sds::Curve::make_rule(keys[i], i));
    }
  }
}

void Slice::relabel_rule(Halfedge_handle h, Sds::Curve c) {
  std::cout << "Relabeling from " << h->curve() << " to " << c << std::endl;
  CGAL_precondition(h->curve().is_inside()==c.is_inside());
  Halfedge_handle hc=h->next()->opposite();
  Sds::Curve ol= h->opposite()->curve();
  h->set_curve(c);
  h->opposite()->set_curve(c.other_side());
  do {
    if (hc->curve() == h->opposite()->curve()) relabel_rule(hc, c.other_side());
    else {
      std::cout << "Not changing " << hc->curve() << std::endl;
    }
    hc= hc->next()->opposite();
  } while (hc != h);
}

Slice::Face_handle Slice::insert_sphere(const T::Sphere_point_3 &ep, 
					T::Key k) {

  Face_handle f;
  Halfedge_handle rvs[4];
  try {
     f= locate_point(ep, k);
     for (unsigned int i=0; i< 4; ++i) {
       rvs[i]= insert_sphere_find_rule_vertex(ep, f, 
					      Sds::Curve::make_rule(k, i));
     }
  } catch (On_edge_exception e) {
  
    if (e.halfedge_handle()->curve().is_rule()) {
      std::cout << "Point hit rule " << e.halfedge_handle()->curve() << std::endl;

      Halfedge_handle h= e.halfedge_handle();
      if (h->curve().is_inside()) h= h->opposite();

      Face_handle pf= h->face();
      Face_handle nf= h->opposite()->face();
      

      int base=0;
      if (h->curve().is_vertical()) ++base;
      int omb=1-base;
      rvs[omb] = insert_sphere_find_rule_vertex(ep, pf, 
						Sds::Curve::make_rule(k, omb));
      rvs[omb+2] = insert_sphere_find_rule_vertex(ep, nf, 
					      Sds::Curve::make_rule(k, omb+2));
      
      //fix
      rvs[base]= h->prev();
      rvs[base+2]= h->opposite()->prev();
      if (!h->curve().is_vertical()) {
	std::swap(rvs[base], rvs[base+2]);
	relabel_rule(h, Sds::Curve(k, Sds::Curve::R_RULE));
	relabel_rule(h->opposite(), Sds::Curve(k, Sds::Curve::L_RULE).other_side());
      } else {
	relabel_rule(h, Sds::Curve(k, Sds::Curve::B_RULE));
	relabel_rule(h->opposite(), Sds::Curve(k, Sds::Curve::T_RULE).other_side());
      }
      f= sds_.remove_rule(h);
      rvs[base]->vertex()->point().replace_rule(Sds::Curve::make_rule(k,base));
      rvs[base+2]->vertex()->point().replace_rule(Sds::Curve::make_rule(k,base+2));
      
    } else {
      std::cerr << "Point hit arc " << e.halfedge_handle()->curve() << std::endl;
      CGAL_assertion(0);
    }
  } catch (On_vertex_exception e) {
    std::cerr << "Point hit vertex " << e.vertex_handle()->point() << std::endl;
    if (e.vertex_handle()->point().type()== Sds::Point::RR) {
      // easy to handle, get 3-4 free ray locations
      // iterate around, check each edge, set the vetex edge if we can
      bool has_rule[4]={false,false, false, false};
      std::vector<Halfedge_handle> hes;
      Halfedge_handle h= e.vertex_handle()->halfedge(), end=h;
      do {
	CGAL_assertion(h->curve().is_rule());
	CGAL_assertion(h->next()->curve().is_rule());
	int base=0;
	if (h->curve().is_vertical()){
	  ++base;
	  if (h->curve().is_inside()) base +=2;
	} else {
	  if (!h->curve().is_inside()) base +=2;
	}
	CGAL_assertion(rvs[base]== Halfedge_handle());
	rvs[base]=h->prev();
	CGAL_assertion(rvs[base]->vertex() == h->opposite()->vertex());
	std::cout << "For " << base << " found vertex " 
		  << rvs[base]->vertex()->point() << std::endl;
	if (h->opposite()->curve().is_inside()) {
	  relabel_rule(h->opposite(), Sds::Curve::make_rule(k, base).other_side());
	} else {
	  relabel_rule(h->opposite(), Sds::Curve::make_rule(k, base));
	}

	rvs[base]->vertex()->point().replace_rule(Sds::Curve::make_rule(k, base));
	std::cout << "Updated to " << rvs[base]->vertex()->point() << std::endl;
	has_rule[base]=true;
	hes.push_back(h);
	h=h->next()->opposite();
      } while (h != end);
      
      for (unsigned int i=0; i< 4; ++i) {
	if (has_rule[i]) continue;
	Halfedge_handle fh=rvs[(i+1)%4];
	CGAL_assertion(rvs[i] == Halfedge_handle());
	CGAL_assertion(has_rule[(i+1)%4]);
	std::cout << "Searching for anchor " << i << " in face ";
	sds_.write_face(fh, std::cout) << std::endl;
	rvs[i]=  insert_sphere_find_rule_vertex(ep, fh->face(), 
						Sds::Curve::make_rule(k,i));
	std::cout << "Got " << rvs[i]->curve() << "--" << rvs[i]->vertex()->point()
		  << "--" << rvs[i]->next()->curve() << std::endl;
      }
      Face_handle f;
      for (unsigned int i=0; i< hes.size(); ++i){
	f=sds_.remove_rule(hes[i]);
      }
    } else {
      CGAL_assertion(0);
    }
  }
  //sds_.audit();

  CGAL_assertion(rvs[3]->face() == rvs[0]->face());
  CGAL_assertion(rvs[0]->face() == rvs[1]->face());
  CGAL_assertion(rvs[1]->face() == rvs[2]->face());
  std::cout << "Creating target..." << std::flush;
  Sds::Halfedge_handle t[4];
  sds_.new_target(k, t);
  std::cout << "done." << std::endl;
  std::cout << "Inserting target..." << std::flush;
  sds_.insert_target(t, rvs);
  std::cout << "done." << std::endl;
  return f;    
}
