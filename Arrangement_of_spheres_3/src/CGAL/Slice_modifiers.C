#include <CGAL/Arrangement_of_spheres_3/Slice.h>


#include <CGAL/IO/qt_debug_examiner_viewer_2.h>



Slice::Halfedge_handle 
Slice::find_rule_vertex(const T::Sphere_point_3 &t, 
			Face_handle f,
			Sds::Curve rule) {
  Halfedge_handle h;
  try {
    h= shoot_rule(t, f, t_.sphere_events(rule.key()).first, rule.rule_direction());
     
      //check_edge_collapse(h->prev());
  } catch (On_vertex_exception e) {
    Vertex_handle v= e.vertex_handle();
    // if it is a rule in the same direction, return it, otherwise
    // pick a random edge
    if (v->point().is_rule_rule() 
	|| v->point().is_sphere_rule() 
	&& v->point().rule_coordinate()== rule.constant_coordinate()) {
      // insert on vertex
      return sds_.find_halfedge(v,f);
    } else {
      // if I am shooting up, make sure I am above the point etc.
      /*Halfedge_handle h0= sds_.find_halfedge(v,f);
      Halfedge_handle h1= h->next();
      if (h0->curve() == h1->curve()) {
	bool cum=false;
	if (rule.is_vertical()) cum = !cum;
	if (h0->curve().arc_index() ==0 || h0->curve().arc_index() ==2) cum= !cum;
	if (cum) {
	  h= h1;
	} else {
	  h= h0;
	}
	
      } else {
	CGAL_assertion(0);
	}*/
      h= sds_.find_halfedge(v,f);
    }
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
    rules[i]= sds_.rule_halfedge(k,Rule_direction(i))->opposite();
    keys[i]= roll_back_rule(ep, rules[i]);
    /*if (keys[i].is_valid()) {
      vs[i]= rules[i]->opposite()->prev();
      CGAL_assertion(vs[i]->vertex() == rules[i]->vertex());
      }*/
    clean_edge(rules[i]);
  }

  //Qt_examiner_viewer_2 *qt= new Qt_examiner_viewer_2();
  *qt_debug_examiner_viewer_2__ << Layer(0);
  draw_rz(qt_debug_examiner_viewer_2__, CGAL::to_double(sim_->current_time()) + .1);
  qt_debug_examiner_viewer_2__->show_everything();
  qt_debug_examiner_viewer_2__->show();
  *qt_debug_examiner_viewer_2__ << std::flush;

  for (unsigned int i=0; i< 4; ++i){
    Halfedge_handle hi= sds_.rule_halfedge(k, Rule_direction(i))->next();
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

  check_merged_faces(sds_.rule_halfedge(k, Rule_direction(0))->face(),
		     sds_.rule_halfedge(k, Rule_direction(2))->face());
  check_merged_faces(sds_.rule_halfedge(k, Rule_direction(1))->face(),
		     sds_.rule_halfedge(k, Rule_direction(3))->face());
  //return Face_handle();
  Halfedge_handle vertices[4];
  Face_handle fr= sds_.remove_target(rules, vertices);
  
  // I guess I don't actually have to extend the rules any more
  // it doesn't  make anything invalid and makes the structure simpler. 
  // Unless there is only halfedge between the two circles
  bool has_split=false;
  for (unsigned int i=0; i< 4; ++i){
    if (keys[i].is_valid()) {
      Rule_direction rd((i+2)%4);
      Sds::Curve c=Sds::Curve::make_rule(keys[i], rd);
      //CGAL_assertion(vertices[i]->vertex()->point().is_sphere_rule());
      Halfedge_handle hv;
      if (!has_split) {
	CGAL_assertion(vertices[rd.index()] != Halfedge_handle());
	std::cout << "using orphaned extremum " 
		  << vertices[i]->vertex()->point() 
		  << " " << vertices[rd.index()]->vertex()->point() << std::endl;
	// I need to handle this
	hv= vertices[rd.index()];
	//vertices[(i+2)%4]=Halfedge_handle();
	keys[rd.index()]= T::Key();
	has_split=true;
      } else {
	std::cout << "Fixing orphaned extremum " 
		  << vertices[i]->vertex()->point() << std::endl;
	// I need to handle this
	// I could just insert it in the new edge since I know it goes there
	hv= find_rule_vertex(ep, vertices[i]->face(),c);
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
						     Vertex_handle vhs[]){
  std::cout << "Point hit rule " << h->curve() << std::endl;
  //Face_handle f;

  if (h->curve().is_inside()) h= h->opposite();
  
  int base=0;
  if (h->curve().is_vertical()) base=3;
  vhs[(base+1)%4] = find_rule_vertex(ep, h->face(), 
				     Sds::Curve::make_rule(k, 
							   Rule_direction((base+1)%4)))->vertex();
  vhs[(base+3)%4] = find_rule_vertex(ep, h->opposite()->face(),
				     Sds::Curve::make_rule(k,
							   Rule_direction((base+3)%4)))->vertex();
  
  
  vhs[base]= h->vertex();
  vhs[(base+2)%4]= h->opposite()->vertex();
  
  clean_edge(h);
  sds_.relabel_rule(h, Sds::Curve::make_rule(k, Rule_direction(base)));
  sds_.relabel_rule(h->opposite(),
		    Sds::Curve::make_rule(k, Rule_direction((base+2)%4)).other_side());
 

  clean_edge(h);
  std::pair<Halfedge_handle, Halfedge_handle> hh= sds_.remove_rule(h);

  return hh.first->face();
}


Slice::Face_handle Slice::insert_sphere_on_rr_prep(const T::Sphere_point_3 &ep, 
						   T::Key k,
						   Vertex_handle v,
						   Vertex_handle vhs[]) {
  // easy to handle, get 3-4 free ray locations

  // halfedge pointing towards v
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
    //CGAL_assertion(rvs[base]== Halfedge_handle());
    hes[base]=h;
    
    // clean up extensions
    clean_edge(h->opposite());
    if (h->opposite()->curve().is_inside()) {
      sds_.relabel_rule(h->opposite(), Sds::Curve::make_rule(k, Rule_direction(base)).other_side());
    } else {
      sds_.relabel_rule(h->opposite(), Sds::Curve::make_rule(k, Rule_direction(base)));
    }
    h=h->next()->opposite();
  } while (h != end);

  // what edge if any is not extending from the vertex
  int missing=-1;

  // delay picking of halfedge until the end in case edges go away

  for (unsigned int i=0; i< 4; ++i) {
    if (hes[i] != Halfedge_handle()){
      vhs[i]= hes[i]->opposite()->vertex();
      continue;
    }
    CGAL_precondition(missing==-1);
    missing=i;
    Halfedge_handle fh=hes[(i+1)%4];
    std::cout << "Searching for anchor " << i << " in face ";
    sds_.write_face(fh, std::cout) << std::endl;
    Halfedge_handle ht=find_rule_vertex(ep, fh->face(),
					Sds::Curve::make_rule(k, Rule_direction(i)));
    vhs[i]= ht->vertex();
    std::cout << "Got " << vhs[i]->point() << std::endl;
  }

  if (missing==-1) {
    missing=0;
    clean_edge(hes[0]);
    sds_.remove_rule(hes[0]);
    hes[0]= Halfedge_handle();
  }
  {
    int curi = (missing+2)%4;
    clean_edge(hes[curi]);
    sds_.remove_rule(hes[curi]);
  }
  clean_edge(hes[(missing+1)%4]);
  clean_edge(hes[(missing+3)%4]);
  Halfedge_handle ht=sds_.remove_redundant_vertex(v->halfedge());
  std::pair<Halfedge_handle, Halfedge_handle> hh = sds_.remove_rule(ht);
  
  Face_handle f= hh.first->face();

  return f;
}




Slice::Face_handle Slice::insert_sphere(const T::Sphere_point_3 &ep, 
					T::Key k) {
  audit();
  Face_handle f;
  // where to put the vertices to attach to
  // the halfedge points to the vertex and is on this face
  Vertex_handle vhs[4];
  try {
    f= locate_point(ep);
  } catch (On_edge_exception e) {
    std::cout << "Point hit edge ";
    sds_.write( e.halfedge_handle(), std::cout) << std::endl;

    if (e.halfedge_handle()->curve().is_rule()) {
      f= insert_sphere_on_rule_prep(ep, k, e.halfedge_handle(), vhs);
    } else {
      Halfedge_handle h= e.halfedge_handle();
      if (!h->curve().is_inside()) h= h->opposite();
      f= h->face();
      std::cout << "Choosing face ";
      sds_.write_face( h, std::cout);
      std::cout << std::endl;
      int start= h->curve().arc_index();
      clean_edge(h);
      vhs[start]=
	sds_.insert_vertex_in_edge(h,
				   Sds::Point(Sds::Curve::make_rule(k, Rule_direction(start)),
					      h->curve()));
      h= sds_.find_halfedge(vhs[start], f)->next();
      vhs[(start+1)%4]=
	sds_.insert_vertex_in_edge(h, 
				   Sds::Point(Sds::Curve::make_rule(k,
								    Rule_direction((start+1)%4)),
					      h->curve()));
    }
  } catch (On_vertex_exception e) {
    std::cerr << "Point hit vertex " << e.vertex_handle()->point() << std::endl;
    if (e.vertex_handle()->point().is_rule_rule()) {
      f= insert_sphere_on_rr_prep(ep, k, e.vertex_handle(), vhs);
    } else if (e.vertex_handle()->point().is_sphere_rule()) {
      // move it to rule
      Halfedge_handle out_rule;
      Halfedge_handle h= e.vertex_handle()->halfedge();
      do {
	if (h->curve().is_rule()){
	  f= insert_sphere_on_rule_prep(ep, k, h, vhs);
	  break;
	}
	h= h->next()->opposite();
      } while (h != e.vertex_handle()->halfedge());
    } else {
      // degeneracy
      CGAL_assertion(0);
    }
  }

  for (unsigned int i=0; i< 4; ++i) {
    if (vhs[i] == Vertex_handle()) {
      vhs[i]= find_rule_vertex(ep, f,
			       Sds::Curve::make_rule(k, Rule_direction(i)))->vertex();
    }
  }
  //sds_.audit();

  /*std::cout << "Creating target..." << std::flush;
  Sds::Halfedge_handle t[4];
  sds_.new_target(k, t);
  std::cout << "done." << std::endl;*/
  Halfedge_handle rvs[4];
  for (unsigned int i=0; i< 4; ++i) {
    rvs[i]= sds_.find_halfedge(vhs[i], f);
  }
  
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
      
  Halfedge_handle hv= find_rule_vertex(ep, f,   ec);
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
					    Halfedge_handle ha,
					    Halfedge_handle hb){
  std::cout << "Calling intersection on ";
  sds_.write(ha, std::cout) << " and ";
  sds_.write(hb, std::cout) << std::endl;
  CGAL_precondition(ha->curve().is_arc());
  CGAL_precondition(hb->curve().is_arc());
  CGAL_precondition(ha->face() == hb->face());
  CGAL_precondition(ha->curve().key() != hb->curve().key());

  // update structure
  clean_edge(ha);
  clean_edge(hb);

  std::pair<Halfedge_handle, Halfedge_handle> nf= sds_.intersect(ha, hb);

  *qt_debug_examiner_viewer_2__ << Layer(0);
  draw_rz(qt_debug_examiner_viewer_2__, CGAL::to_double(sim_->current_time()) + .1);
  qt_debug_examiner_viewer_2__->show_everything();
  qt_debug_examiner_viewer_2__->show();
  *qt_debug_examiner_viewer_2__ << std::flush;

  sds_.audit();

  std::cout << "inside face is ";
  sds_.write_face(nf.first, std::cout) << std::endl;
  std::cout << "first outside face is ";
  sds_.write_face(nf.first->opposite(), std::cout) << std::endl;
  std::cout << "second outside face is ";
  sds_.write_face(nf.second->opposite(), std::cout) << std::endl;
  std::cout << "first corner face is ";
  sds_.write_face(nf.first->opposite()->prev()->opposite(), std::cout) << std::endl;
  std::cout << "second corner face is ";
  sds_.write_face(nf.second->opposite()->prev()->opposite(), std::cout) << std::endl;

  // update edges
  check_edge_collapse(nf.first);
  check_edge_collapse(nf.first->opposite()->prev());
  check_edge_collapse(nf.first->next()->opposite()->next());
  check_edge_collapse(nf.second);
  check_edge_collapse(nf.second->opposite()->prev());
  check_edge_collapse(nf.second->next()->opposite()->next());
  check_edge_face(nf.first->opposite());
  check_edge_face(nf.second->opposite());
  audit();
  return nf.first->face();
}

Slice::Face_handle Slice::unintersect_spheres(const T::Event_point_3 &ep,
					      Halfedge_handle ha,
					      Halfedge_handle hb) {
  // find face
  check_merged_faces(ha->opposite()->prev()->opposite()->face(),
		     hb->opposite()->prev()->opposite()->face());
  clean_edge(ha->opposite()->prev());
  clean_edge(ha->opposite()->prev()->opposite()->prev());
  clean_edge(hb->opposite()->prev());
  clean_edge(hb->opposite()->prev()->opposite()->prev());
  std::pair<Halfedge_handle, Halfedge_handle> rp= sds_.unintersect(ha->face());
  check_edge_collapse(rp.first);
  check_edge_collapse(rp.second);
  audit();
  return rp.first->face();
}



Slice::Face_handle Slice::collapse_edge(const T::Event_point_3 &ep,
					Halfedge_handle rule,
					Halfedge_handle c,
					Halfedge_handle ne){
  CGAL_precondition(sds_.degree(rule->vertex()) ==3);
  CGAL_precondition(c->curve().is_arc());

  roll_back_rule(ep, rule->opposite());

  {
    *qt_debug_examiner_viewer_2__ << Layer(0);
    draw_rz(qt_debug_examiner_viewer_2__, CGAL::to_double(sim_->current_time()) + .1);
    qt_debug_examiner_viewer_2__->show_everything();
    qt_debug_examiner_viewer_2__->show();
    *qt_debug_examiner_viewer_2__ << std::flush;
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

  /*{
    Qt_examiner_viewer_2 *qt= new Qt_examiner_viewer_2();
    draw_rz(qt, CGAL::to_double(sim_->current_time()) + .1);
    qt->show_everything();
    qt->show();
    }*/

  if (sds_.is_redundant(vo)) {
    clean_edge(vo->halfedge());
    clean_edge(vo->halfedge()->next());
    Halfedge_handle h= sds_.remove_redundant_vertex(vo->halfedge());
    check_edge_collapse(h);
  }

  check_edge_collapse(sds_.next_edge_on_curve(ne->opposite()));
  return Face_handle();
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
  Halfedge_handle hv= find_rule_vertex(ep, of, rule->curve());
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
