#include <CGAL/Arrangement_of_spheres_3/Irrational_cross_section_insertion.h> 

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Irrational_cross_section_insertion CGAL_AOS3_TARG ::CS::Face_handle 
Irrational_cross_section_insertion CGAL_AOS3_TARG::insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
							       CGAL_AOS3_TYPENAME CS::Face_handle f) {

  // where to put the vertices to attach to
  // the halfedge points to the vertex and is on this face
  CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4];

  return finish_insert(k, f, vhs);
}



CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Irrational_cross_section_insertion CGAL_AOS3_TARG ::CS::Face_handle 
Irrational_cross_section_insertion CGAL_AOS3_TARG::finish_insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
								      CGAL_AOS3_TYPENAME CS::Face_handle f,
								      CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4]) {
  std::cout << "Finishing insert into face ";
  cs_.write(f, std::cout);
  std::cout << "\n with vertices ";
  for (unsigned int i=0; i< 4; ++i) {
    if (vhs[i] != CGAL_AOS3_TYPENAME CS::Vertex_handle()){
      cs_.write(vhs[i], std::cout) << " ";
    } else {
      std::cout << "null ";
    }
  }
  std::cout << std::endl;
  //P::cs_.audit(true);
  ICSR icsr(tr_, cs_);
  for (unsigned int i=0; i< 4; ++i) {
    if (vhs[i] == CGAL_AOS3_TYPENAME CS::Vertex_handle()) {
      vhs[i]= icsr.find_rule_vertex(tr_.sphere_events(k).first, f,
				    CS::Curve::make_rule(k, Rule_direction(i)))->vertex();
    }
  }
  //P::cs_.audit(true);
  CGAL_AOS3_TYPENAME CS::Halfedge_handle rvs[4];
  for (unsigned int i=0; i< 4; ++i) {
    rvs[i]= cs_.find_halfedge(vhs[i], f);
  }
  
  std::cout << "Inserting target..." << std::flush;
  CGAL_AOS3_TYPENAME CS::Halfedge_handle cvs[4];
  cs_.new_target(k, cvs);
    
  //std::vector<CGAL_AOS3_TYPENAME CS::Halfedge_handle> out;
  const CGAL_AOS3_TYPENAME CS::Curve curves[4]={ CS::Curve::make_rule(k, Rule_direction(0)),
						 CS::Curve::make_rule(k, Rule_direction(1)).opposite(),
						 CS::Curve::make_rule(k, Rule_direction(2)).opposite(),
						 CS::Curve::make_rule(k, Rule_direction(3))};
  cs_.star_face(rvs, rvs+4, curves, curves+4);
  CGAL_AOS3_TYPENAME CS::Halfedge_handle inward_edges[4];
  for (unsigned int i=0; i< 4; ++i) {
    inward_edges[i]= rvs[i]->next();
  }
  cs_.expand_vertex(cvs, cvs+4, inward_edges, inward_edges+4 );
  //cs_.stitch_in(cvs, cvs+4, rvs, curves);
  std::cout << "done." << std::endl;
  for (unsigned int i=0; i< 4; ++i){
    //P::cse_.check_edge_collapse(rvs[i]->next());
    //if (P::cs_.event(rvs[i]) == CGAL_AOS3_TYPENAME CS::Event_key()) {
    //P::cse_.check_edge_collapse(rvs[i]);
    //}
    // if (P::cs_.event(rvs[i]->next()->opposite()->next()) == CGAL_AOS3_TYPENAME CS::Event_key()) {
    //P::cse_.check_edge_collapse(rvs[i]->next()->opposite()->next());
    //}
    //P::cse_.check_edge_face(rvs[i]->next()->next());
    //check_edge_collapse(rvs[i]->next()->next());
  }
  // clean up if there are degeneracies
  for (unsigned int i=0; i< 4; ++i) {
    if (vhs[i]->point().is_rule_rule()) {
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h= vhs[i]->halfedge();
      if (h->opposite()->vertex()->point().is_rule_rule() && h->opposite()->vertex()->point().is_finite()) {
	if (cs_.next_edge_on_rule(cs_.cross_edge(h)->opposite()) != CGAL_AOS3_TYPENAME CS::Halfedge_handle()) {
	  std::cout << "edge ";
	  cs_.write(h, std::cout) << " is superfluous" <<  std::endl;
	  cs_.join_face(h, true);
	  --i; // just restart the loop
	  continue;
	}
      }
    }
  }

  CGAL_LOG_WRITE(Log::LOTS, cs_.write(LOG_STREAM));
  cs_.audit();
  return cs_.a_halfedge(k)->face();    
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Irrational_cross_section_insertion CGAL_AOS3_TARG ::CS::Face_handle 
Irrational_cross_section_insertion CGAL_AOS3_TARG::insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
							       CGAL_AOS3_TYPENAME CS::Halfedge_handle h) {

  std::cout << "Point hit edge ";
  cs_.write( h, std::cout) << std::endl;
  CGAL_AOS3_TYPENAME CS::Face_handle f;
  CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4];
  ICSR icsr(tr_, cs_);
  if (h->curve().is_rule()) {
    icsr.roll_back_rule(tr_.sphere_events(k).first, h);
    icsr.roll_back_rule(tr_.sphere_events(k).first, h->opposite());
    int hi=h->curve().halfedge_direction().index();
    int hoi=h->opposite()->curve().halfedge_direction().index();
    CGAL_assertion((hi+2)%4==hoi);
    vhs[h->curve().halfedge_direction().index()] = h->vertex();
    vhs[h->opposite()->curve().halfedge_direction().index()]= h->opposite()->vertex();
    int hp1=(h->curve().halfedge_direction().index()+1) %4;
    vhs[hp1]= icsr.find_rule_vertex(tr_.sphere_events(k).first, h->face(),
				    CS::Curve::make_rule(k, Rule_direction(hp1)))->vertex();
    int hop1=(h->curve().halfedge_direction().index()+3) %4;
    vhs[hop1]= icsr.find_rule_vertex(tr_.sphere_events(k).first, h->opposite()->face(),
				     CS::Curve::make_rule(k, Rule_direction(hop1)))->vertex();
    cs_.audit(true);
    f= cs_.join_face(h, false);
  } else {
    //CGAL_AOS3_TYPENAME CS::Halfedge_handle h= e;
    if (!h->curve().is_inside()) h= h->opposite();
    f= h->face();
    std::cout << "Choosing face ";
    cs_.write( h, std::cout);
    std::cout << std::endl;
    int start= h->curve().arc_index();
    //P::cse_.clean_edge(h);
    CGAL_AOS3_TYPENAME CS::Vertex_handle nvh=cs_.new_vertex(CGAL_AOS3_TYPENAME CS::Point( CS::Curve::make_rule(k,Rule_direction(start)),h->curve()));
    h = cs_.insert_vertex(nvh, h);
    vhs[start] = h->vertex();
    //h= cs_.find_halfedge(vhs[start], f)->next();
    h=h->next();

    CGAL_AOS3_TYPENAME CS::Vertex_handle nvh2= cs_.new_vertex(CGAL_AOS3_TYPENAME CS::Point(CS::Curve::make_rule(k,
												Rule_direction((start+1)%4)),
												h->curve()));
    vhs[(start+1)%4]=
      cs_.insert_vertex(nvh, h)->vertex();
  }
  return finish_insert(k,f,vhs);
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Irrational_cross_section_insertion CGAL_AOS3_TARG ::CS::Face_handle 
Irrational_cross_section_insertion CGAL_AOS3_TARG ::insert(CGAL_AOS3_TYPENAME Traits::Sphere_3_key k,
							       CGAL_AOS3_TYPENAME CS::Vertex_handle v) {
  std::cerr << "Point hit vertex " << v->point() << std::endl;
  CGAL_AOS3_TYPENAME CS::Vertex_handle vhs[4];
  CGAL_AOS3_TYPENAME CS::Face_handle f;
  ICSR icsr(tr_, cs_);
  if (v->point().is_rule_rule()) {
      
    std::vector<CGAL_AOS3_TYPENAME CS::Halfedge_handle> to_remove;
    CGAL_AOS3_TYPENAME CS::Halfedge_handle h=v->halfedge()->opposite();
    CGAL_AOS3_TYPENAME CS::Face_handle missing_face;
    do {
      icsr.roll_back_rule(tr_.sphere_events(k).first, h);
      to_remove.push_back(h);
      //P::cse_.clean_edge(h);
      CGAL_AOS3_TYPENAME CS::Rule_direction rd= h->curve().halfedge_direction();
      std::cout << "Curve " << h->curve() << " has direction " << rd.to_str() 
		<< " and index " << rd.index() << std::endl;
      vhs[rd.index()]= h->vertex();
      if (h->prev()->curve().halfedge_direction() == h->curve().halfedge_direction()) {
	missing_face= h->face();
      }
      h= h->prev()->opposite();
    } while (h != v->halfedge()->opposite()); 
      
    for (unsigned int i=0; i< 4; ++i) {
      if (vhs[i] == CGAL_AOS3_TYPENAME CS::Vertex_handle()) {
	vhs[i]= icsr.find_rule_vertex(tr_.sphere_events(k).first, missing_face,
				      CS::Curve::make_rule(k, Rule_direction(i)))->vertex();

      }

    }
    cs_.audit(true);
    CGAL_AOS3_TYPENAME CS::Vertex_handle iv= to_remove.front()->opposite()->vertex();
    CGAL_AOS3_TYPENAME CS::Face_handle fh= cs_.unstar_face(iv, false);

      /*cs_.snip_out(to_remove.begin(),
					     to_remove.end(), false);*/
    //cs_.delete_component(iv);


    /*for (unsigned int i=0; i< 4; ++i) {
      if (vhs[i] == CGAL_AOS3_TYPENAME CS::Vertex_handle()) {
      vhs[i]= find_rule_vertex(tr_.sphere_events(k).first, f,
      CGAL_AOS3_TYPENAME CS::Curve::make_rule(k, Rule_direction(i)))->vertex();
      }
      }*/
    return finish_insert(k, fh, vhs);
  } else if (v->point().is_sphere_rule()) {
    // move it to rule
    CGAL_AOS3_TYPENAME CS::Halfedge_handle out_rule;
    CGAL_AOS3_TYPENAME CS::Halfedge_handle h= v->halfedge();
    do {
      if (h->curve().is_rule()){
	return insert(k,h);
      }
      h= h->next()->opposite();
    } while (h != v->halfedge());
    CGAL_assertion(0);
    return CGAL_AOS3_TYPENAME CS::Face_handle();
  } else {
    // degeneracy
    return insert(k, v->halfedge()->face());
  }
}

CGAL_AOS3_END_INTERNAL_NAMESPACE
