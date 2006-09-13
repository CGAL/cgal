#include <CGAL/Arrangement_of_spheres_3/Slice.h>


#include <CGAL/IO/Qt_examiner_viewer_2.h>

/*Slice::Face_handle Slice::insert_sphere(const T::Sphere_point_3 &cp) {
  T::Key k=t_.new_sphere(cp.sphere());
  halfedges_.resize(halfedges_.size()+1);
  boost::array<Vertex_handle, 4> vf={{0}};
  rule_certificates_.push_back(vf);
  return insert_sphere(t_.sphere_events(k).first, k);
  }*/



Slice::Halfedge_handle 
Slice::find_rule_vertex(const T::Sphere_point_3 &ep, 
			Face_handle f,
			Sds::Curve rule) {
  Halfedge_handle h;
  try {
      h= shoot_rule(ep, f, rule);
     
      //check_edge_collapse(h->prev());
  } catch (On_vertex_exception e) {
    Vertex_handle v= e.vertex_handle();
    // pick random edge incident on face and insert it there
    h= sds_.find_halfedge(v,f);

    /*if (v->point().is_sphere_extremum()) {
      bool inside=false;
      {
	Halfedge_handle h= v->halfedge();
	do {
	  if (h->curve().is_arc()) {
	    if (h->face()==f && h->curve().is_inside()
		|| h->opposite()->face() == f && !h->curve().is_inside()) {
	      inside=true;
	      break;
	    }
	  }
	  h=h->next()->opposite();
	} while (h != v->halfedge());
      }
     
      if (inside) {
	// v=v;
      } else {
	Halfedge_handle h= v->halfedge();
	do {
	  h=h->next()->opposite();
	} while (!h->curve().is_rule());
	clean_edge(h);
	v=sds_.insert_vertex_in_edge(h, Sds::Point(rule, h->curve()));
	//check_edge_collapse(h);
	CGAL_assertion(h->prev()->vertex()==v);
	//check_edge_collapse(h->prev());
      } 
    } else if (v->point().is_rule_rule()) {
      
    } else if (v->point().is_sphere_rule())*/
  }
  clean_edge(h);
  Vertex_handle v=sds_.insert_vertex_in_edge(h, Sds::Point(rule, h->curve()));
  //check_edge_collapse(h);
  CGAL_assertion(h->prev()->vertex()==v);
  return sds_.find_halfedge(v,f);
}


Slice::T::Key Slice::roll_back_rule(const T::Sphere_point_3 &ep,
				    Halfedge_handle cur) {
  std::cout << "Rolling back, rolling back..." << cur->curve()
	    << " at " << cur->vertex()->point() << std::endl;
    //<< " with last of ";
  /*if ( last==Halfedge_handle()) std::cout << "null";
    else std::cout << last->curve();*/
  T::Key k;
  Halfedge_handle next= sds_.next_edge_on_curve(cur);
  if (next != Halfedge_handle()) {
    if (next->curve().key() == cur->curve().key()) {
      std::cout << "Recursing..." << std::endl;
      k=roll_back_rule(ep, next);
    } else {
      std::cout << "terminating on..." << next->curve() << std::endl;
      k= next->curve().key();
    }
  } else {
    if (cur->vertex()->point().is_sphere_rule()) {
      // we hit another sphere;
      /*CGAL_assertion(cur->vertex()->point().sphere(0).key() 
	!= cur->vertex()->point().rule(0).key());*/
      if (cur->vertex()->point().is_sphere_extremum()
	  && !cur->next()->curve().is_inside()) {
	k= cur->next()->curve().key();
      } 
    } 
  }
   

  if (!cur->opposite()->vertex()->point().is_sphere_extremum()) {
    if (k != T::Key() && cur->curve().key() != k) {
      cur->curve().flip_rule(k);
    } else {
      Halfedge_handle vhh= cur->opposite()->prev();
      //sds_.write(vhh, std::cout);
      CGAL_assertion(vhh->vertex() == cur->vertex());
      rotate_rule(ep, cur->opposite());
      /*std::cout << "Checking ";
      sds_.write(vhh, std::cout);
      std::cout << std::endl;
      if (sds_.is_redundant(vhh->vertex())) {
	clean_edge(vhh);
	clean_edge(vhh->next());
	Halfedge_handle h= sds_.remove_redundant_vertex(vhh);
	check_edge_collapse(h);
      } else {
	std::cout << "Vertex " << vhh->vertex()->point() 
		  << " does not need to be removed." << std::endl;
		  }*/
    }
  }
  return k;
}							   

Slice::Face_handle Slice::erase_sphere(const T::Sphere_point_3 &ep,
				       T::Key k) {
  
  audit();
  //T::Key k= f->halfedge()->curve().key();

  Face_handle f;

  f= sds_.a_halfedge(k)->face();
 
  // check correctness of f
  Halfedge_handle rules[4];
  Halfedge_handle h=f->halfedge();
  int deg=0;
  do {
    ++deg;
    CGAL_assertion(h->curve().key()==k);
    CGAL_assertion(h->curve().is_arc());
    CGAL_assertion(h->curve().is_inside());
    /*Halfedge_handle rule= h->opposite()->next();
      int index= rule->curve().rule_index();*/
    /*CGAL_assertion(rules[index]
		   == Halfedge_handle());
		   rules[index] =rule;*/
    h= h->next();
  } while (h != f->halfedge());
  CGAL_assertion(deg==4);
  // roll in each until I have a target in a face
  
  T::Key keys[4];
  //Halfedge_handle vs[4];
  for (unsigned int i=0; i< 4; ++i){
    rules[i]= sds_.rule_halfedge(k,i)->opposite();
    keys[i]= roll_back_rule(ep, rules[i]);
    /*if (keys[i].is_valid()) {
      vs[i]= rules[i]->opposite()->prev();
      CGAL_assertion(vs[i]->vertex() == rules[i]->vertex());
      }*/
  }

  Qt_examiner_viewer_2 *qt= new Qt_examiner_viewer_2();
  draw_rz(qt, CGAL::to_double(sim_->current_time()) + .1);
  qt->show_everything();
  qt->show();

  for (unsigned int i=0; i< 4; ++i){
    Halfedge_handle hi= sds_.rule_halfedge(k, i)->next();
    if (sds_.event(hi) != Event_key()) {
      std::cerr << "ERROR " << sds_.event(hi) 
		<< " on edge ";
      sds_.write(hi, std::cerr);
      std::cerr << std::endl;
    }
    //CGAL_assertion(sds_.event(hi) == Event_key());
    /*if (sds_.event(hi->opposite()->prev())
	!= sim_->null_event()) {
      std::cerr << "ERROR " << sds_.event(hi->opposite()->prev())
		<< " on edge ";
      sds_.write(hi->opposite()->prev(), std::cerr);
      std::cerr << std::endl;
      }*/
    clean_edge(hi->opposite()->prev());
  }

  check_merged_faces(sds_.rule_halfedge(k, 0)->face(),
		     sds_.rule_halfedge(k, 2)->face());
  check_merged_faces(sds_.rule_halfedge(k, 1)->face(),
		     sds_.rule_halfedge(k, 3)->face());
  //return Face_handle();
  Halfedge_handle vertices[4];
  Face_handle fr= sds_.remove_target(rules, vertices);
  
  // I guess I don't actually have to extend the rules any more
  // it doesn't  make anything invalid and makes the structure simpler. 
  // Unless there is only halfedge between the two circles
  bool has_split=false;
  for (unsigned int i=0; i< 4; ++i){
    if (keys[i].is_valid()) {
      Sds::Curve c=Sds::Curve::make_rule(keys[i], (i+2)%4);
      //CGAL_assertion(vertices[i]->vertex()->point().is_sphere_rule());
      Halfedge_handle hv;
      if (!has_split) {
	CGAL_assertion(vertices[(i+2)%4] != Halfedge_handle());
	std::cout << "using orphaned extremum " 
		  << vertices[i]->vertex()->point() 
		  << " " << vertices[(i+2)%4]->vertex()->point() << std::endl;
	// I need to handle this
	hv= vertices[(i+2)%4];
	//vertices[(i+2)%4]=Halfedge_handle();
	keys[(i+2)%4]= T::Key();
	has_split=true;
      } else {
	std::cout << "Fixing orphaned extremum " 
		  << vertices[i]->vertex()->point() << std::endl;
	// I need to handle this
	// I could just insert it in the new edge since I know it goes there
	hv= find_rule_vertex(ep, vertices[i]->face(), 
			     c);
	//vertices[i]->vertex()->point().replace_rule(c);
      }
      Halfedge_handle nh;
      if (i==0 || i == 3) {
	nh=sds_.split_face(hv,vertices[i], c);
      } else {
	nh=sds_.split_face(vertices[i], hv, c);
      }
      check_edge_collapse(nh);
    } 
  }
  for (unsigned int i=0; i< 4; ++i){ 
    if (sds_.is_redundant(vertices[i]->vertex())) {
      clean_edge(vertices[i]);
      clean_edge(vertices[i]->next());
      Halfedge_handle h= sds_.remove_redundant_vertex(vertices[i]);
      if (sds_.is_redundant(h->vertex())
	  || sds_.is_redundant(h->opposite()->vertex())) {
	// skipping since it is adjacent to another vertex which will be removed
	// seems like there should be a better way
	std::cout << "Skipping ephemeral edge ";
	sds_.write(h, std::cout);
	std::cout << std::endl;
      } else {
	check_edge_collapse(h);
      }
    }
  }

  check_reduced_face(fr);

  //audit();
  return fr;
}

void Slice::relabel_rule(Halfedge_handle h, Sds::Curve c) {
  std::cout << "Relabeling from " << h->curve() << " to " << c << std::endl;
  clean_edge(h);
  CGAL_precondition(h->curve().is_inside()==c.is_inside());
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

Slice::Halfedge_handle Slice::check_remove_redundant(Halfedge_handle h) {
  if (sds_.is_redundant(h->vertex())) {
    CGAL_assertion(h->curve()
		   == h->next()->curve());
    // fails for degeneracies, I'll fix that elsewhere (in the degeneracy handler)
    clean_edge(h);
    clean_edge(h->next());
    Halfedge_handle hn= sds_.remove_redundant_vertex(h);
    check_edge_collapse(hn);
    return hn;
  } else {
    std::cout << "Not removing vertex " <<h->vertex()->point() << std::endl;
    return h;
  }
}

Slice::Face_handle Slice::insert_sphere_on_rule_prep(const T::Sphere_point_3 &ep, 
						     T::Key k,
						     Halfedge_handle h,
						     Halfedge_handle rvs[]){
  std::cout << "Point hit rule " << h->curve() << std::endl;
  Face_handle f;

  if (h->curve().is_inside()) h= h->opposite();
  
  Face_handle pf= h->face();
  Face_handle nf= h->opposite()->face();
  
  
  int base=0;
  if (h->curve().is_vertical()) ++base;
  int omb=1-base;
  rvs[omb] = find_rule_vertex(ep, pf, 
			      Sds::Curve::make_rule(k, omb));
  rvs[omb+2] = find_rule_vertex(ep, nf, 
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
  clean_edge(h);
  std::pair<Halfedge_handle, Halfedge_handle> hh= sds_.remove_rule(h);
  check_remove_redundant(hh.first);
  check_remove_redundant(hh.second);
  //f= sds_.remove_rule(h);
  //rvs[base]->vertex()->point().replace_rule(Sds::Curve::make_rule(k,base));
  //rvs[base+2]->vertex()->point().replace_rule(Sds::Curve::make_rule(k,base+2));

  

  return hh.first->face();
}


Slice::Face_handle Slice::insert_sphere_on_rr_prep(const T::Sphere_point_3 &ep, 
						   T::Key k,
						   Vertex_handle v,
						   Halfedge_handle rvs[]) {
  // easy to handle, get 3-4 free ray locations
  // iterate around, check each edge, set the vetex edge if we can
  bool has_rule[4]={false,false, false, false};
  Halfedge_handle hes[4];
  Halfedge_handle h= v->halfedge(), end=h;
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
    
    //rvs[base]->vertex()->point().replace_rule(Sds::Curve::make_rule(k, base));
    std::cout << "Updated to " << rvs[base]->vertex()->point() << std::endl;
    has_rule[base]=true;
    hes[base]=h;
    h=h->next()->opposite();
  } while (h != end);
      
  for (unsigned int i=0; i< 4; ++i) {
    if (has_rule[i]) continue;
    Halfedge_handle fh=rvs[(i+1)%4];
    CGAL_assertion(rvs[i] == Halfedge_handle());
    CGAL_assertion(has_rule[(i+1)%4]);
    std::cout << "Searching for anchor " << i << " in face ";
    sds_.write_face(fh, std::cout) << std::endl;
    rvs[i]=  find_rule_vertex(ep, fh->face(), 
					    Sds::Curve::make_rule(k,i));
    std::cout << "Got " << rvs[i]->curve() << "--" << rvs[i]->vertex()->point()
	      << "--" << rvs[i]->next()->curve() << std::endl;
  }
  Face_handle f;
  for (unsigned int i=0; i< 4; ++i){
    if (has_rule[i]) {
      Sds::Curve r= Sds::Curve::make_rule(k, i);
      if (i==1 || i ==2) r=r.other_side();
      relabel_rule(hes[i]->opposite(), r);
      clean_edge(hes[i]);
      std::pair<Halfedge_handle, Halfedge_handle> hh=sds_.remove_rule(hes[i]);
      check_remove_redundant(hh.first);
      check_remove_redundant(hh.second);
      f= hh.first->face();
    }
  }
  return f;
}

Slice::Face_handle Slice::insert_sphere_on_ss_prep(const T::Sphere_point_3 &ep, 
						   T::Key k,
						   Vertex_handle h,
						   Halfedge_handle rvs[]){
  CGAL_assertion(0);
  return Face_handle();
}


Slice::Face_handle Slice::insert_sphere_on_arc_prep(const T::Sphere_point_3 &ep, 
						    T::Key k,
						    Halfedge_handle h,
						    Halfedge_handle rvs[]){
  std::cerr << "Point hit arc " << h->curve() << std::endl;
  Face_handle f;
  if (h->curve().is_inside()) h=h->opposite();
  //CGAL::Comparison_result cmp= ep.compare(t_.center(e.halfedge_handle()->curve().key()), T::Coordinate_index(2));
  
  Vertex_handle v0, v1;
  v0= sds_.insert_vertex_in_edge(h, Sds::Point::make_special(k));
  v1= sds_.insert_vertex_in_edge(h, Sds::Point::make_special(k));
  Halfedge_handle h0= h->prev()->prev();
  CGAL_assertion(h0->vertex()==v0);
  CGAL_assertion(h0->vertex() == v0);
  Halfedge_handle h1= h->prev();
  CGAL_assertion(h1->vertex()==v1);
  int start=2;
  int step=1;
  //if ( cmp!=CGAL::SMALLER) {
  f= h->face();
  //start=2;
  //step=1;
  /*} else {
    start=1;
    f= h->opposite()->face();
    step=3;
    h0= h0->opposite()->prev();
    h1= h1->opposite()->prev();
    }*/

  start+= h->curve().arc_index();
  std::cout << "Arc index is " << h->curve().arc_index() << std::endl;
  rvs[start%4]= h0;
  rvs[(start+step)%4]= h1;
  h0->vertex()->point() = Sds::Point(h->curve(),
				     Sds::Curve::make_rule(k, start));

  h1->vertex()->point() = Sds::Point(h->curve(),
				     Sds::Curve::make_rule(k, (start+step)%4));
  std::cout << "v0 is " << h0->vertex()->point() << std::endl;
  std::cout << "v1 is " << h1->vertex()->point() << std::endl;
  for (unsigned int i=0; i< 4; ++i) {
    if (rvs[i] == Halfedge_handle()) {
      rvs[i]=find_rule_vertex(ep, f, 
					    Sds::Curve::make_rule(k, i));
    }
  }
  sds_.write_face(h0, std::cout) << std::endl;
  return f;
} 


Slice::Face_handle Slice::insert_sphere(const T::Sphere_point_3 &ep, 
					T::Key k) {
  audit();
  Face_handle f;
  Halfedge_handle rvs[4];
  try {
     f= locate_point(ep, k);
     for (unsigned int i=0; i< 4; ++i) {
       rvs[i]= find_rule_vertex(ep, f, 
					      Sds::Curve::make_rule(k, i));
     }
  } catch (On_edge_exception e) {
  
    if (e.halfedge_handle()->curve().is_rule()) {
      f= insert_sphere_on_rule_prep(ep, k, e.halfedge_handle(), rvs);
    } else {
      f= insert_sphere_on_rule_prep(ep, k, e.halfedge_handle(), rvs);
    }
  } catch (On_vertex_exception e) {
    std::cerr << "Point hit vertex " << e.vertex_handle()->point() << std::endl;
    if (e.vertex_handle()->point().is_rule_rule()) {
      f= insert_sphere_on_rr_prep(ep, k, e.vertex_handle(), rvs);
    } else if (e.vertex_handle()->point().is_sphere_rule()) {
      // move it outside. if the rule is outside put it on it
      // otherwise insert it on the arc
      Halfedge_handle out_rule;
      Halfedge_handle h= e.vertex_handle()->halfedge();
      do {
	if (h->curve().is_rule() && h->next()->curve().is_arc() 
	    &&  !h->next()->curve().is_inside()) {
	  out_rule=h;
	  break;
	}
	h= h->next()->opposite();
      } while (h != e.vertex_handle()->halfedge());

      if (out_rule != Halfedge_handle()) {
	f= insert_sphere_on_rule_prep(ep,k, h, rvs);
      } else {
	f= insert_sphere_on_arc_prep(ep, k, h, rvs);
      }
    } else {
      f= insert_sphere_on_ss_prep(ep, k, e.vertex_handle(), rvs);
    }
  }
  //sds_.audit();

  CGAL_assertion(rvs[3]->face() == rvs[0]->face());
  CGAL_assertion(rvs[0]->face() == rvs[1]->face());
  CGAL_assertion(rvs[1]->face() == rvs[2]->face());
  /*std::cout << "Creating target..." << std::flush;
  Sds::Halfedge_handle t[4];
  sds_.new_target(k, t);
  std::cout << "done." << std::endl;*/
  std::cout << "Inserting target..." << std::flush;
  sds_.insert_target(k, rvs);
  std::cout << "done." << std::endl;
  for (unsigned int i=0; i< 4; ++i){
    check_edge_collapse(rvs[i]->next());
    if (sds_.event(rvs[i]) == Event_key()) {
      check_edge_collapse(rvs[i]);
    }
    if (sds_.event(rvs[i]->next()->opposite()->next()) == Event_key()) {
      check_edge_collapse(rvs[i]->next()->opposite()->next());
    }
    check_edge_face(rvs[i]->next()->next());
    //check_edge_collapse(rvs[i]->next()->next());
  }
  
  audit();
  return f;    
}


// makes rule go away
Slice::Halfedge_handle Slice::rotate_rule(const T::Event_point_3 &ep,
					  Halfedge_handle rule) {
  CGAL_precondition(rule->curve().is_rule());
  CGAL_precondition(rule->vertex()->point().is_rule_rule());
  
  std::cout << "Rotating rule " << rule->curve() << " about vertex "
	    << rule->vertex()->point() << std::endl;

  Halfedge_handle oe;
  Face_handle f;
  Sds::Curve ec;
  if (rule->next()->curve().key() == rule->curve().key()) {
    oe= rule;
    f= rule->face();
    ec= rule->opposite()->prev()->curve();
  } else {
    oe= rule->opposite()->prev();
    f= rule->opposite()->face();
    ec= rule->next()->opposite()->curve();
  }
      
  Halfedge_handle hv= find_rule_vertex(ep, f,  ec);
  std::cout << "New edge supported by " << ec 
	    << " from " << oe->vertex()->point() << " to " 
	    << hv->vertex()->point() << std::endl;
  Halfedge_handle nh=sds_.split_face(oe, hv, ec);
  check_edge_collapse(nh);
  
  Vertex_handle t= rule->vertex();
  std::cout << "Removing rule ";
  sds_.write(rule, std::cout);
  std::cout << std::endl;
  clean_edge(rule);
    
  Halfedge_handle th= rule->opposite()->prev();
  CGAL_assertion(th->vertex() == t);

  check_merged_faces(rule->face(), rule->opposite()->face());
  std::pair<Halfedge_handle, Halfedge_handle>hp= sds_.remove_rule(rule);
  check_remove_redundant(hp.first);
  check_remove_redundant(hp.second);
  return hv;
}



Slice::Face_handle Slice::intersect_spheres(const T::Event_point_3 &t,
					    T::Key k, T::Key l){
  // find face(s)
  typedef std::pair<Halfedge_handle, Halfedge_handle> EP;
  EP ep;
  {
    std::vector<EP> shared_faces;
    Halfedge_handle kh= sds_.a_halfedge(k);
    Halfedge_handle khs=kh;
    do {
      Halfedge_handle oh= kh->opposite()->face()->halfedge();
      do {
	if (oh->curve().key() == l && oh->curve().is_arc()) {
	  shared_faces.push_back(EP(kh->opposite(), oh));
	}
	oh= oh->next();
      } while (oh != kh->opposite()->face()->halfedge());
      
      kh= sds_.next_edge_on_curve(kh);
      CGAL_assertion(kh != Halfedge_handle());
    } while (kh != khs);
    // find edges
    CGAL_assertion(!shared_faces.empty());
    if (shared_faces.size() !=1) {
      CGAL_assertion(0);
    } 
    ep= shared_faces[0];
  }

  // update structure
  clean_edge(ep.first);
  clean_edge(ep.second);
  Face_handle nf= sds_.intersect(ep.first, ep.second);

  Qt_examiner_viewer_2 *qt= new Qt_examiner_viewer_2();
  draw_rz(qt, CGAL::to_double(sim_->current_time()) + .1);
  qt->show_everything();
  qt->show();

  sds_.audit();

  // update edges
  if (true) {
    Halfedge_handle ch= nf->halfedge();
    for (int i=0; i< 2;++i) {
      check_edge_collapse(ch->opposite()->next());
      check_edge_collapse(ch->opposite()->prev());
      check_edge_face(ch->opposite());
      ch= ch->next();
      check_reduced_face(ch->opposite()->prev()->opposite()->face());
    }
  }
  audit();
  return nf;
}

Slice::Face_handle Slice::unintersect_spheres(const T::Event_point_3 &ep,
					      T::Key k, T::Key l){
  // find face
  Face_handle f;
  {
    Halfedge_handle kh= sds_.a_halfedge(k);
    Halfedge_handle khs=kh;
    do {
      if (kh->next()->curve().key() == l 
	  && kh->next()->next() == kh){
	f= kh->face();
      } else if (kh->opposite()->next()->curve().key() ==l
	  && kh->opposite()->next()->next() == kh->opposite()) {
	f= kh->opposite()->face();
      }
      kh= sds_.next_edge_on_curve(kh);
      CGAL_assertion(kh != Halfedge_handle());
    } while (kh != khs);
    // find edges
    //CGAL_assertion(!shared_faces.empty());
  }
  // clear edge certs
  Halfedge_handle h0= f->halfedge();
  Halfedge_handle h1= h0->next();

  clean_edge(h0->opposite()->prev());
  clean_edge(h0->opposite()->next());
  clean_edge(h1->opposite()->prev());
  clean_edge(h1->opposite()->next());

  // modify structure
  std::pair<Halfedge_handle, Halfedge_handle> nes= sds_.unintersect(f);
  // create new edge certs
  check_edge_collapse(nes.first);
  check_edge_collapse(nes.second);

  Qt_examiner_viewer_2 *qt= new Qt_examiner_viewer_2();
  draw_rz(qt, CGAL::to_double(sim_->current_time()) + .1);
  qt->show_everything();
  qt->show();

  audit();

  return nes.first->face();
}



Slice::Face_handle Slice::collapse_edge(const T::Event_point_3 &ep,
					Halfedge_handle rule,
					Halfedge_handle c,
					Halfedge_handle ne){
  CGAL_precondition(sds_.degree(rule->vertex()) ==3);
  CGAL_precondition(c->curve().is_arc());

  roll_back_rule(ep, rule->opposite());

  {
    Qt_examiner_viewer_2 *qt= new Qt_examiner_viewer_2();
    draw_rz(qt, CGAL::to_double(sim_->current_time()) + .1);
    qt->show_everything();
    qt->show();
  }

  check_merged_faces(rule->face(), rule->opposite()->face());

  clean_edge(c);
  Vertex_handle vc= sds_.insert_vertex_in_edge(c, Sds::Point(c->curve(),
							    rule->curve()));
  Vertex_handle vo= rule->opposite()->vertex();

  clean_edge(rule->next());
  clean_edge(rule->opposite()->prev());

  Halfedge_handle cp= c->prev();

  /*rule=*/ sds_.move_rule(rule, ne, cp); 
  check_edge_collapse(cp);

  {
    Qt_examiner_viewer_2 *qt= new Qt_examiner_viewer_2();
    draw_rz(qt, CGAL::to_double(sim_->current_time()) + .1);
    qt->show_everything();
    qt->show();
  }

  if (sds_.is_redundant(vo)) {
    clean_edge(vo->halfedge());
    clean_edge(vo->halfedge()->next());
    Halfedge_handle h= sds_.remove_redundant_vertex(vo->halfedge());
    check_edge_collapse(h);
  }

  check_edge_collapse(sds_.next_edge_on_curve(ne->opposite()));
}

Slice::Halfedge_handle Slice::uncollapse_edge(const T::Event_point_3 &ep,
					      Halfedge_handle rule,
					      Halfedge_handle c,
					      Halfedge_handle ne) {
  CGAL_assertion(sds_.degree(rule->vertex()) ==3);
  Halfedge_handle xe;
  Halfedge_handle base=c; CGAL_assertion(0);
  if (base->face() == rule->face()) {
    xe= base->next();
  } else {
    xe= base->prev();
  }
  Face_handle of= xe->opposite()->face();
  Halfedge_handle hv= find_rule_vertex(ep, of,  rule->curve());
  check_edge_collapse(hv);
  check_edge_collapse(hv->next());
  Halfedge_handle ov= rule->opposite()->prev();
  Halfedge_handle ob= rule->next();
  if (ob== base) {
    ob= rule->opposite()->prev();
  }
  sds_.move_rule(rule->opposite(), base, hv); 
  if (sds_.is_redundant(ov->vertex())) {
    clean_edge(ov);
    clean_edge(ov->next());
    CGAL_assertion(0);
    //Halfedge_handle h= sds_.remove_redundant_vertex(ov->vertex());
    //check_edge_collapse(h);
  }
  check_edge_collapse(ob);
  CGAL_assertion(0);
  return ov;
}
