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
  Halfedge_handle v= insert_vertex(Sds::Point(rule, h->curve()), h);
  //check_edge_collapse(h);
  //CGAL_assertion(h->prev()->vertex()==v);
  CGAL_assertion(v->face() == f);
  return v;
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
	nh=sds_.split_face(c, hv,vertices[i]);
      } else {
	nh=sds_.split_face(c, vertices[i], hv);
      }
      check_edge_collapse(nh);
    } 
  }
  for (unsigned int i=0; i< 4; ++i){ 
    if (check_remove_vertex(vertices[i]) != vertices[i]) {
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
  Halfedge_handle hn= h->next();
  sds_.merge_faces(h);

  return hn->face();
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
      sds_.relabel_rule(h->opposite(),
			Sds::Curve::make_rule(k,
					      Rule_direction(base)).other_side());
    } else {
      sds_.relabel_rule(h->opposite(), 
			Sds::Curve::make_rule(k, Rule_direction(base)));
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
					Sds::Curve::make_rule(k,
							      Rule_direction(i)));
    vhs[i]= ht->vertex();
    std::cout << "Got " << vhs[i]->point() << std::endl;
  }

  if (missing==-1) {
    missing=0;
    clean_edge(hes[0]);
    sds_.merge_faces(hes[0]);
    hes[0]= Halfedge_handle();
  }
  {
    int curi = (missing+2)%4;
    clean_edge(hes[curi]);
    sds_.merge_faces(hes[curi]);
  }
  clean_edge(hes[(missing+1)%4]);
  clean_edge(hes[(missing+3)%4]);
  Halfedge_handle ht=check_remove_vertex(v->halfedge());
  Halfedge_handle hn= ht->next();
  sds_.merge_faces(ht);
  
  
  return hn->face();
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
      h =
	sds_.insert_vertex(Sds::Point(Sds::Curve::make_rule(k,
							    Rule_direction(start)),
				      h->curve()), h);
      vhs[start] = h->vertex();
      //h= sds_.find_halfedge(vhs[start], f)->next();
      h=h->next();
      vhs[(start+1)%4]=
	insert_vertex(Sds::Point(Sds::Curve::make_rule(k,
						       Rule_direction((start+1)%4)),
		      h->curve()), h)->vertex();
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
  
  sds_.audit();

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


/*
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
							     rule->curve()))->vertex();
  Vertex_handle vo= rule->opposite()->vertex();

  clean_edge(rule->next());
  clean_edge(rule->opposite()->prev());

  Halfedge_handle cp= c->prev();

  sds_.move_rule(rule, ne, cp); 
  check_edge_collapse(cp);

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
*/



Slice::Halfedge_handle Slice::collapse_rule(const T::Event_point_3 &ep,
					    Halfedge_handle rule,
					    Halfedge_handle circle1) {
  CGAL_precondition(rule->vertex()->point().is_sphere_extremum());
  CGAL_precondition(!circle1->curve().is_inside());

  roll_back_rule(ep, rule->opposite());
  //CGAL_assertion(ok== T::Key());
  //bool right=false;
  Halfedge_handle arc0, arc1, circle0, circle2, c;
  if (rule->vertex() != circle1->vertex()) {
    // RHS case
    arc1= circle1->next();
    circle2= arc1->opposite()->next();
    circle0= rule->opposite()->prev();
    arc0=circle1->opposite()->prev();
    c= insert_vertex(circle1->vertex()->point(),
		 circle0);
  } else {
    //right=true;
    arc1=circle1->prev();
    circle2= arc1->opposite()->prev();
    circle0= rule->next();
    arc0=circle1->opposite()->next();

    c= insert_vertex(circle2->vertex()->point(),
				      circle0);
  }
  circle0= Halfedge_handle();//
  Halfedge_handle r0= insert_vertex(Sds::Point(rule->curve(),//
					       arc1->curve()),
				    rule->opposite());
  Halfedge_handle r1= r0->next();//
  rule= Halfedge_handle();//
  
  arc0 = move_edge_target(arc0, c->opposite()->prev());//
  split_face(arc0->curve(), c, r0);//
  move_edge_target(arc1, r0);//
  if (!r1->vertex()->point().is_sphere_extremum()) {//
    Halfedge_handle v= r1->opposite()->prev();
    merge_faces(r1);
    check_remove_vertex(v);
  }
  
  check_remove_vertex(circle1);
 
  return r0;
}

Slice::Halfedge_handle Slice::uncollapse_rule(const T::Event_point_3 &ep,
					      Halfedge_handle r,
					      Halfedge_handle c1) {
  Halfedge_handle c2, c0, a0, r0, r1, a2, a1;
  Vertex_handle vh;
  if (c1->vertex() == r->vertex()) {
    c2= r->next();
    a2= c1->prev();
    a1= c1->opposite()->next();
    c0= a1->opposite()->prev();
    a0= c1->opposite()->next()->opposite();
  } else {
    c1= r->opposite()->next();
    c2= r->prev();
    a2= r->next()->opposite();
    a1= r->opposite()->prev();
    c0= a1->opposite()->next();
    a0= c1->opposite()->prev();
  }
  // handle degeneracy later
  if (sds_.degree(r->opposite()->vertex()) ==3) {
    Halfedge_handle v= find_rule_vertex(ep, c0->face(), r->curve());
    split_face(r->curve(), r->prev()->opposite()->prev(), v);
  }
  Halfedge_handle c= insert_vertex(a1->vertex()->point(), c2);
  move_edge_target(a0, c);
  merge_faces(a1);
  check_remove_vertex(vh->halfedge());
  move_edge_target(a2, c);
  check_remove_vertex(r->opposite());
}
