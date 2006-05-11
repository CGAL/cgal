#ifndef CGAL_NEF3_REFLEX_EDGE_SEARCHER_H
#define CGAL_NEF3_REFLEX_EDGE_SEARCHER_H

#include<CGAL/Nef_3/SNC_decorator.h>

CGAL_BEGIN_NAMESPACE

template<typename Nef_>
class Reflex_edge_searcher : public Modifier_base<typename Nef_::SNC_structure> {

  typedef Nef_                                            Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_structure          SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>              SNC_decorator;  

  typedef typename SNC_structure::Vertex_handle           Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle         Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle        Halffacet_handle;
  typedef typename SNC_structure::SHalfedge_handle        SHalfedge_handle;
  typedef typename SNC_structure::SHalfloop_handle        SHalfloop_handle;
  typedef typename SNC_structure::SFace_handle            SFace_handle;

  typedef typename SNC_structure::Volume_iterator         Volume_iterator;
  typedef typename SNC_structure::Shell_entry_iterator    Shell_entry_iterator;

  typedef typename SNC_structure::Vector_3                Vector_3;
  typedef typename SNC_structure::Point_3                 Point_3;

  typedef typename SNC_structure::Sphere_segment          Sphere_segment;

  typedef typename std::deque<Halfedge_handle>             Edge_list;
 public:
  typedef typename std::deque<Halfedge_handle>::iterator   Reflex_edge_iterator;
  typedef Edge_list                                        Container;

 private:
  struct Reflex_edge_visitor {
    
    Edge_list& reflex_edges;
    
    Reflex_edge_visitor(Edge_list& e) : reflex_edges(e) {}
    
    void visit(Vertex_handle v) const {}
    void visit(Halfedge_handle e) const {}
    void visit(Halffacet_handle f) const {}
    void visit(SHalfloop_handle sl) const {}
    void visit(SFace_handle sf) const {}
    void visit(SHalfedge_handle se) const {
      Halfedge_handle e = se->source();
      //      std::cerr << "is reflex edge?" << std::endl;
      //      std::cerr << "  e " << e->source()->point() 
      //		<< "->" << e->twin()->source()->point() << std::endl;
      //      std::cerr << "  marks " << se->incident_sface()->mark() << ", "
      //		<< se->incident_sface()->volume()->mark() << std::endl;
      if(e->source()->point() > e->twin()->source()->point())
	return;
      SHalfedge_handle se2 = se->sprev()->twin();
      CGAL_assertion(se->source() == se2->source());

      if(se2 == se) {
	//	std::cerr << "problem found " << std::endl;
	reflex_edges.push_back(e);
	return;
      }
      //      std::cerr << " se1 " << se->circle() << std::endl;
      //      std::cerr << " se2 " << se2->circle() << std::endl;
      Vector_3 vec1 = e->point() - CGAL::ORIGIN;
      Vector_3 vec2 = se->circle().orthogonal_vector();
      Point_3 p_comp = CGAL::ORIGIN + cross_product(vec1,vec2);
      CGAL_assertion(Sphere_segment(e->point(),
				    e->twin()->point(),
				    se->twin()->circle()).has_on(p_comp));
      //      std::cerr << "test side of " << p_comp 
      //		<< " in relation to " << se->circle() 
      //		<< ": " << se2->circle().oriented_side(p_comp) << std::endl;
      if(se2->circle().oriented_side(p_comp) == ON_NEGATIVE_SIDE)
	reflex_edges.push_back(e);


    }
  };
  
 public:
  bool first;
  Volume_iterator c;
  Volume_iterator c_end;
  Edge_list reflex_edges;

  Reflex_edge_searcher() : first(true) {}

  void operator()(SNC_structure& snc) {
    reflex_edges.clear();
    
    Reflex_edge_visitor rev(reflex_edges);
    SNC_decorator D(snc);
    CGAL_forall_volumes(c, snc) {
      //      std::cerr << "new volume " << std::endl;
      if(c->mark())
	for(Shell_entry_iterator shi=c->shells_begin(); shi!=c->shells_end(); ++shi) {
	  //	  std::cerr << "new shell " << std::endl;
	  D.visit_shell_objects(SFace_handle(shi),rev);
	}
    }
  }

  void append(Halfedge_handle e) { 
    if(e->twin()->source()->point() < e->source()->point())
      e = e->twin();
    reflex_edges.push_back(e);
  }
  
  Edge_list& get_container() { return reflex_edges; }
  Reflex_edge_iterator reflex_edges_begin() { return reflex_edges.begin(); }
  Reflex_edge_iterator reflex_edges_end() { return reflex_edges.end(); }
};

CGAL_END_NAMESPACE
#endif // CGAL_NEF3_REFLEX_EDGE_SEARCHER_H
