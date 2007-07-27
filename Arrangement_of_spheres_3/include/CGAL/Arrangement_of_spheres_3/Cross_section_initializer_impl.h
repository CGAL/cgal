#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_initializer.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_arrangement.h>


CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE


CGAL_AOS3_TEMPLATE 
Cross_section_initializer CGAL_AOS3_TARG::Cross_section_initializer(Combinatorial_cross_section CGAL_AOS3_TARG &cs,
								    const Traits &tr): cs_(cs),
										       traits_(tr){
}

CGAL_AOS3_TEMPLATE 
void Cross_section_initializer CGAL_AOS3_TARG::operator()(CGAL_AOS3_TYPENAME Traits::FT z ) {
  cs_.clear();
  typedef Cross_section_arrangement CGAL_AOS3_TARG Arr;
  Arr arr(traits_.sphere_3s_begin(),
				traits_.sphere_3s_end(),z,1000000);

  
  for (CGAL_AOS3_TYPENAME Arr::Face_iterator fit= arr.faces_begin(); 
       fit != arr.faces_end(); ++fit){
    new_face(fit->begin(), fit->end());
  }
  /*if (gen_certs) {
    initialize_certificates();
  }
    
  sds_.set_is_building(false);*/

  /*if (gen_certs) {
    for (Sds::Halfedge_iterator it= sds_.halfedges_begin(); 
	 it != sds_.halfedges_end(); ++it){
      if (it->curve().is_inside() && it->curve().is_finite()) check_edge_collapse(it);
    }
    for (Sds::Face_iterator it= sds_.faces_begin(); it != sds_.faces_end(); ++it){
      Halfedge_handle c= it->halfedge();
      do {
	check_edge_face(c);
	c= c->next();
      } while (c != it->halfedge());
    }
  }
  if (gen_certs) {
    audit();
  }*/

  for (CGAL_AOS3_TYPENAME std::map<Edge, CGAL_AOS3_TYPENAME CS::Halfedge_handle>::iterator it= unmatched_hedges_.begin();
       it != unmatched_hedges_.end(); ++it){
    //std::cout << "Searching for next for ";
    //write(it->second, std::cout) << std::endl;
    it->second->set_face(cs_.infinite_face());
    cs_.infinite_face()->set_halfedge(it->second);
    CGAL_AOS3_TYPENAME CS:: Halfedge_handle c=it->second->opposite()->prev()->opposite();
    CGAL_AOS3_TYPENAME CS:: Vertex_handle v= it->second->vertex();
    while (c->prev() != CGAL_AOS3_TYPENAME CS::Halfedge_handle()){
      //write( c, std::cout) << std::endl;
      CGAL_AOS3_TYPENAME CS::Vertex_handle vo= c->opposite()->vertex();
      CGAL_assertion(v==vo);
      c= c->prev()->opposite();
    }
    //std::cout << "Found ";
    //write(c, std::cout) << std::endl;
    c->set_prev(it->second);
    it->second->set_next(c);
  }
  unmatched_hedges_.clear();
  
  for (CGAL_AOS3_TYPENAME CS::Halfedge_iterator it= cs_.halfedges_begin(); 
       it != cs_.halfedges_end(); ++it){
    if (it->curve().is_arc() && it->curve().key().is_input()
	&& it->curve().is_inside()) {
      CGAL_AOS3_TYPENAME CS::Halfedge_handle fit= cs_.next_edge_on_curve(it);
      if (fit->curve() != it->curve()) {
	int ai= (it->curve().arc_index()+1)%4;
	cs_.halfedges(it->curve().key())[ai]=it;
      }
    }
  }
  
}


CGAL_AOS3_TEMPLATE
Cross_section_initializer CGAL_AOS3_TARG::~Cross_section_initializer() {

}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Vertex_handle
Cross_section_initializer CGAL_AOS3_TARG::new_vertex_cached(CGAL_AOS3_TYPENAME CS::Point p) {
  CGAL_precondition(p.is_valid());
  //std::cout << "Creating point " << p << std::endl;
  CGAL_AOS3_TYPENAME CS::HDS::Vertex v;
  v.point()=p;
  points_[p]=cs_.hds_.vertices_push_back(v);
  return points_[p];
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Combinatorial_cross_section CGAL_AOS3_TARG::Halfedge_handle 
Cross_section_initializer CGAL_AOS3_TARG::new_halfedge(CGAL_AOS3_TYPENAME CS::Point s, 
						       CGAL_AOS3_TYPENAME CS::Curve ff, 
						       CGAL_AOS3_TYPENAME CS::Point f) {
  //std::cout << "Creating edge " << s << " -- " << ff << " -- " << f << std::endl;
  CGAL_precondition(points_.find(s) != points_.end());
  if (points_.find(f) == points_.end()) {
    new_vertex_cached(f);
  }

  CGAL_AOS3_TYPENAME CS::Halfedge_handle h;
  Edge ep(s, ff, f);
  if (unmatched_hedges_.find(ep) != unmatched_hedges_.end()){
    h= unmatched_hedges_[ep];
    unmatched_hedges_.erase(ep);
    CGAL_assertion(h->is_border());
    //std::cout << "matched" << std::endl;
  } else {
    //std::cout << "unmatched" << std::endl;
    h= cs_.new_halfedge(ff);
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

CGAL_AOS3_END_INTERNAL_NAMESPACE

