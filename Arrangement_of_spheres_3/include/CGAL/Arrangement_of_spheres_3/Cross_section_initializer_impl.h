#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_initializer.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_arrangement.h>

#include <CGAL/IO/Qt_examiner_viewer_2.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE


CGAL_AOS3_TEMPLATE 
Cross_section_initializer CGAL_AOS3_TARG::Cross_section_initializer(Combinatorial_cross_section CGAL_AOS3_TARG &cs,
								    const Traits &tr): cs_(cs),
										       traits_(tr){
  //CGAL_precondition(cs_.number_of_vertices() ==4);
  //CGAL_precondition(cs_.number_of_faces() == 2);
}


CGAL_AOS3_TEMPLATE 
void
Cross_section_initializer CGAL_AOS3_TARG::copy_in(const Arr &arr,
						  std::vector<Halfedge_handle> &ties) {
  for (CGAL_AOS3_TYPENAME Arr::Face_iterator fit= arr.faces_begin(); 
       fit != arr.faces_end(); ++fit){
    new_face(fit->begin(), fit->end(), arr);
  }
  ++cs_.num_components_;
  
  std::vector<CGAL_AOS3_TYPENAME Arr::Halfedge_const_handle> outer_edges;
  arr.outer_face(outer_edges);
  Halfedge_handle h= new_face(outer_edges.begin(), outer_edges.end(), arr);
  Halfedge_handle c=h;
  do {
    if (cs_.degree(c->vertex()) ==1) {
      ties.push_back(c);
    }
    
    c=c->prev();
  } while (c != h);
  cs_.visit_component(h);
}


CGAL_AOS3_TEMPLATE 
void Cross_section_initializer CGAL_AOS3_TARG::operator()(CGAL_AOS3_TYPENAME Traits::FT z ) {
  cs_.clear();
  cs_.num_components_=1;
  cs_.set_number_of_spheres(traits_.number_of_sphere_3s());
  Arr arr(traits_,z);

  cs_.set_has_boundary(true);
  
  std::vector<CGAL_AOS3_TYPENAME Arr::Halfedge_const_handle> outer_edges;
  arr.outer_points(outer_edges);

  if (!outer_edges.empty()) {
    std::vector<Curve > points[4];
    int offset=-1;
    Curve cc, cn;
    do {
      ++offset;
      cc= arr.curve(outer_edges[offset]);
      cn= arr.curve(outer_edges[(offset+1)%outer_edges.size()]);
      CGAL_LOG(Log::SOME, "Curves are " << cc << " and " << cn << std::endl);
    }  while (!(cc.rule_outward_direction() == Rule_direction(0) 
		&& cn.rule_outward_direction() == Rule_direction(3)));
    ++offset;
		
    for (unsigned int i=0; i< outer_edges.size(); ++i) {
      CGAL_AOS3_TYPENAME Arr::Halfedge_const_handle hc=outer_edges[(offset+i)%outer_edges.size()];
      //PP hpt= ;
      Curve c= arr.curve(hc);
      CGAL_assertion(hc->target()->degree());
      //Sphere_3_key x= pt.rule_key(plane_coordinate(0));
      //Sphere_3_key y= pt.rule_key(plane_coordinate(1));
      CGAL_LOG(Log::SOME, "Putting " << c << " into bin " 
	       << c.rule_outward_direction() << std::endl);
      points[c.rule_outward_direction().index()].push_back(c);
    }
    
        
    Halfedge_handle starth= cs_.infinite_face()->halfedge();
    CGAL_LOG(Log::SOME, "Starth is ");
    CGAL_LOG_WRITE(Log::SOME, cs_.write(starth, LOG_STREAM) << std::endl);

    while (starth->curve().rule_outward_direction() != Rule_direction(3)) {
      starth=starth->next();
    }
    
    starth=starth->opposite();
    CGAL_LOG_WRITE(Log::SOME, cs_.write(starth->face(), LOG_STREAM) << std::endl);
    std::vector<Halfedge_handle> handles;
    std::vector<Curve> curves;
 
    for (int i=0; i< 4; ++i ) {
      std::reverse(points[i].begin(), points[i].end()); 
      for (unsigned int j=0; j < points[i].size(); ++j) {
	CGAL_LOG(Log::SOME, "Curve is " << points[i][j] << " and Starth is ");
	CGAL_LOG_WRITE(Log::SOME, cs_.write(starth, LOG_STREAM) << std::endl);
	
	Point npt(points[i][j],
		  starth->curve());
	CGAL_LOG(Log::SOME, "Inserting " << npt << " into edge ");
	CGAL_LOG_WRITE(Log::SOME, cs_.write(starth, LOG_STREAM) << std::endl);
	Vertex_handle vh= cs_.new_vertex(npt);
	cs_.insert_vertex(vh, starth);
	handles.push_back(starth->prev());
	//starth=starth->next();
	curves.push_back(points[i][j]);
	CGAL_LOG(Log::SOME, "Curve/handle pair is " << curves.back() << ": ");
	CGAL_LOG_WRITE(Log::SOME, cs_.write(handles.back(), LOG_STREAM) << std::endl);
      }
      starth=starth->next();
    }
    
    CGAL_LOG_WRITE(Log::SOME, cs_.write(starth->face(), LOG_STREAM) 
		   << " is inside face" << std::endl);
    
    Vertex_handle vh=cs_.star_face(handles.begin(), handles.end(), curves.begin(), curves.end());

    std::vector<Halfedge_handle> inner_halfedges;
    copy_in(arr, inner_halfedges);

    std::vector<Halfedge_handle> outer_halfedges;
    {
      int offset=0;
      while (handles[offset]->next()->curve()
	     != inner_halfedges[0]->opposite()->curve()) {
	++offset;
      }
      for (unsigned int i=0; i< handles.size(); ++i) {
	outer_halfedges.push_back(handles[(i+offset)%handles.size()]->next());
      }
    }

    cs_.expand_vertex(inner_halfedges.begin(), inner_halfedges.end(),
		      outer_halfedges.begin(), outer_halfedges.end());


    for (CGAL_AOS3_TYPENAME CS::Halfedge_iterator it= cs_.halfedges_begin(); 
	 it != cs_.halfedges_end(); ++it){
      if (it->curve().is_arc() && it->curve().key().is_input()
	  && it->curve().is_inside()) {
	CGAL_AOS3_TYPENAME CS::Halfedge_handle fit= cs_.next_edge_on_circle(it);
	if (fit->curve() != it->curve()) {
	  //int ai= (it->curve().arc_index()+1)%4;
	  cs_.halfedges_[it->curve().key().input_index()]=it;
	}
      }
      if (it->curve().is_inside()) {
	cs_.v_.on_new_edge(it);
      }
    }


  }
 
  CGAL_LOG_WRITE(Log::LOTS,cs_.write(LOG_STREAM));
  
  cs_.audit();


}


CGAL_AOS3_TEMPLATE
Cross_section_initializer CGAL_AOS3_TARG::~Cross_section_initializer() {
 
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Vertex_handle
Cross_section_initializer CGAL_AOS3_TARG::find_vertex(CGAL_AOS3_TYPENAME Arr::Vertex_const_handle hc,
						      const Arr &arr) {
  //CGAL_precondition(p.is_valid());
  //std::cout << "Creating point " << p << std::endl;
  Vertex_handle vh;
  if (points_.find(hc) == points_.end()) {
    vh=cs_.new_vertex(arr.point(hc));
    points_[hc]=vh;
  } 
  return points_[hc];
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle 
Cross_section_initializer CGAL_AOS3_TARG::find_halfedge(CGAL_AOS3_TYPENAME Arr::Halfedge_const_handle hc,
						       const Arr &arr) {
  //std::cout << "Creating edge " << s << " -- " << ff << " -- " << f << std::endl;
  //CGAL_precondition(points_.find(hc->source()) != points_.end());
  
  CGAL_AOS3_TYPENAME CS::Halfedge_handle h;
  
  if (unmatched_hedges_.find(hc) != unmatched_hedges_.end()){
    h= unmatched_hedges_[hc];
    //unmatched_hedges_.erase(hc);
    CGAL_assertion(h->is_border());
    CGAL_LOG_WRITE(Log::SOME, cs_.write(h, LOG_STREAM) << " was found"
		   << std::endl);
    //std::cout << "matched" << std::endl;
  } else {
    //std::cout << "unmatched" << std::endl;
    h= cs_.new_halfedge(arr.curve(hc));
    //h->set_inside(inside);
    //h->opposite()->set_inside(!inside);
    unmatched_hedges_[hc->twin()]= h->opposite();
    unmatched_hedges_[hc]= h;
    find_vertex(hc->target(), arr)->set_halfedge(h);
    find_vertex(hc->source(), arr)->set_halfedge(h->opposite());
    h->set_vertex( find_vertex(hc->target(), arr));
    h->opposite()->set_vertex( find_vertex(hc->source(), arr));
    CGAL_LOG_WRITE(Log::SOME, cs_.write(h, LOG_STREAM) << " is new halfedge"
		   << std::endl);
    CGAL_LOG_WRITE(Log::SOME, cs_.write(h->opposite(), LOG_STREAM) << " is new unmatched"
		   << std::endl);
  }
  return h;
}

CGAL_AOS3_END_INTERNAL_NAMESPACE

