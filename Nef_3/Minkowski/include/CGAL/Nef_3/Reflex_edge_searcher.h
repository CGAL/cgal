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
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator
                                  SHalfedge_around_svertex_circulator;

  typedef typename SNC_structure::Vector_3                Vector_3;
  typedef typename SNC_structure::Point_3                 Point_3;

  typedef typename SNC_structure::Sphere_point            Sphere_point;
  typedef typename SNC_structure::Sphere_circle           Sphere_circle;
  typedef typename SNC_structure::Sphere_segment          Sphere_segment;

  typedef typename std::deque<SHalfedge_handle>           SEdge_list;

 public:
  typedef typename std::deque<SHalfedge_handle>::iterator   Reflex_sedge_iterator;
  typedef SEdge_list                                        Container;

 private:
  struct Reflex_edge_visitor {
    
    SEdge_list& pos;
    SEdge_list& neg;
    Sphere_point dir;
      
    Reflex_edge_visitor(SEdge_list& p, SEdge_list& n, Sphere_point dir_in) 
        : pos(p), neg(n), dir(dir_in) {}
    
    void visit(Vertex_handle v) const {}
    void visit(Halfedge_handle e) const {}
    void visit(Halffacet_handle f) const {}
    void visit(SHalfloop_handle sl) const {}
    void visit(SFace_handle sf) const {}
    void visit(SHalfedge_handle se) const {
      int isrse = Reflex_edge_searcher::is_reflex_sedge(se, dir);
//      std::cerr << "isrse " << isrse << std::endl;
      if((isrse&1)==1) pos.push_back(se);
      if((isrse&2)==2) neg.push_back(se);
      if((isrse&2)==2) {
          se->source()->mark()=false;
          se->source()->twin()->mark()=false;
      }
    }
  };
  
 public:
  SEdge_list pos;
  SEdge_list neg;
  Sphere_point dir;

  Reflex_edge_searcher(Sphere_point dir_in) 
      : dir(dir_in) {}

/*
  bool is_reflex_edge(Halfedge_handle e) {
    SHalfedge_around_svertex_circulator svc(e->out_sedge()),
        end(svc);
    CGAL_For_all(svc,end)
      if(svc->incident_sface()->mark() &&
        is_reflex_sedge(svc))
          return true;
    return false;
  }
*/
  int is_reflex_sedge(SHalfedge_handle se) {
    return is_reflex_sedge(se, dir);
  }

  static int is_reflex_sedge(SHalfedge_handle se, Sphere_point dir) {

    Halfedge_handle e = se->source();
//    std::cerr << "is reflex edge?" << std::endl;
//    std::cerr << "  e " << e->source()->point() 
//              << "->" << e->twin()->source()->point() << std::endl;
//    std::cerr << "  marks " << se->incident_sface()->mark() << std::endl;
    if(e->source()->point() > e->twin()->source()->point())
	return 0;
    if(e->point() == dir || e->twin()->point() == dir)
        return 0;

    SHalfedge_handle se2 = se->sprev()->twin();
    
    if(se2 == se)
	return 3;

//    std::cerr << " se1 " << se->circle() << std::endl;
//    std::cerr << " se2 " << se2->circle() << std::endl;
    Vector_3 vec1 = e->point() - CGAL::ORIGIN;
    Vector_3 vec2 = se->circle().orthogonal_vector();
    Sphere_point sp1 = CGAL::ORIGIN + cross_product(vec2,vec1);
    if(se2->circle().oriented_side(sp1) != ON_POSITIVE_SIDE)
      return 0;
    
    int result = 0;
    Sphere_circle cp(e->point(), dir);
//    std::cerr << " cp " << se2->circle() << std::endl;
    Vector_3 vec3 = cp.orthogonal_vector();
    Sphere_point sp3 = CGAL::ORIGIN + cross_product(vec3,vec1);
    CGAL::Oriented_side os1 = se->circle().oriented_side(sp3);
    CGAL::Oriented_side os2 = se2->circle().oriented_side(sp3);

    if(os1 == ON_POSITIVE_SIDE || os2 == ON_NEGATIVE_SIDE)
        result |= 1;
    if(os1 == ON_NEGATIVE_SIDE || os2 == ON_POSITIVE_SIDE)
        result |= 2;

//    std::cerr << "result " << result << std::endl;

/*
    vec2 = se2->circle().orthogonal_vector();
    Sphere_point sp2 = CGAL::ORIGIN + cross_product(vec2,vec1);
    Sphere_segment s(sp1, sp2, Sphere_circle(sp2,sp1));

    CGAL_assertion(s.is_long());
    


    int result = 0;
    Sphere_point sp(dir);
    if(s.has_on(sp) && s.source() != sp && s.target() != sp)
      result |= 1;
    
    sp = sp.antipode();
    if(s.has_on(sp) && s.source() != sp && s.target() != sp)
      result |= 2;
*/
//    CGAL_assertion(result != 0);
    return result;
  }

  void operator()(SNC_structure& snc) {
    pos.clear();
    neg.clear();
    
    Reflex_edge_visitor rev(pos,neg,dir);
    SNC_decorator D(snc);
    Volume_iterator c;
    CGAL_forall_volumes(c, snc) {
      //      std::cerr << "new volume " << std::endl;
      if(c->mark())
	for(Shell_entry_iterator shi=c->shells_begin(); shi!=c->shells_end(); ++shi) {
	  //	  std::cerr << "new shell " << std::endl;
	  D.visit_shell_objects(SFace_handle(shi),rev);
	}
    }
  }

  void handle_new_edge(Halfedge_handle e) { 
    if(e->twin()->source()->point() < e->source()->point())
      e = e->twin();
    SHalfedge_around_svertex_circulator svc(e->out_sedge()), send(svc);
    CGAL_For_all(svc, send) {
      int isrse = is_reflex_sedge(svc, dir);
      if(isrse == 0) continue;
      if((isrse&1==1)) pos.push_back(svc);
      if((isrse&2==2)) neg.push_back(svc);
      break;
    }
  }

  SEdge_list& get_positive_rsedges() { return pos; }
  SEdge_list& get_negative_rsedges() { return neg; }
  Reflex_sedge_iterator positive_rsedges_begin() { return pos.begin(); }
  Reflex_sedge_iterator positive_rsedges_end() { return pos.end(); }
  Reflex_sedge_iterator negative_rsedges_begin() { return neg.begin(); }
  Reflex_sedge_iterator negative_rsedges_end() { return neg.end(); }
};

CGAL_END_NAMESPACE
#endif // CGAL_NEF3_REFLEX_EDGE_SEARCHER_H
