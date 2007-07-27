#ifndef CGAL_ARRANGEMENT_OF_SPHERES_CROSS_SECTION_INIT_H
#define CGAL_ARRANGEMENT_OF_SPHERES_CROSS_SECTION_INIT_H

#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE
CGAL_AOS3_TEMPLATE
struct Cross_section_initializer {
  CGAL_AOS3_TRAITS;
public:
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
  Cross_section_initializer (CS& cs, const Traits& tr);

  void operator()(CGAL_AOS3_TYPENAME Traits::FT z);

  struct Edge {
    
    Edge(CGAL_AOS3_TYPENAME CS::Point s, 
	 CGAL_AOS3_TYPENAME CS::Curve su,
	 CGAL_AOS3_TYPENAME CS::Point t): sup_(su), s_(s), t_(t){}
    
    bool operator<(const Edge &o) const {
      if (sup_ < o.sup_) return true;
      else if (o.sup_ < sup_) return false;
      else if (s_ < o.s_) return true;
      else if (o.s_ < s_) return false;
      else return t_ < o.t_;
    }
    
    CGAL_AOS3_TYPENAME CS::Curve sup_;
    CGAL_AOS3_TYPENAME CS::Point s_, t_;
  };
  

  template <class It>
  void new_face(It b, It e) {
    // vt must be a pair, ick 
    CGAL_assertion(std::distance(b,e)>=2);
    It em1=e;
    --em1;
    CGAL_AOS3_TYPENAME CS::Point lp= em1->second;
    if (points_.find(lp)==points_.end()) new_vertex_cached(lp);
      
    CGAL_AOS3_TYPENAME CS::Face_handle f= cs_.hds_.faces_push_back(CS::HDS::Face());
    It c=b;
    CGAL_AOS3_TYPENAME CS::Halfedge_handle le, fst;
    while (c != e) {
      CGAL_AOS3_TYPENAME CS::Point cp= c->second;
      CGAL_AOS3_TYPENAME CS::Halfedge_handle h= new_halfedge(lp, c->first, cp);
      if (c==b) fst=h;
      h->set_face(f);
      f->set_halfedge(h);
      lp=cp;
      if (le != CS::Halfedge_handle()) {
	le->set_next(h);
	h->set_prev(le);
	CGAL_assertion(le->next()->opposite()->vertex() == le->vertex());
      }
      ++c;
      le=h;
    }
    le->set_next(fst);
    fst->set_prev(le);
  }

  ~Cross_section_initializer();
private:

  CGAL_AOS3_TYPENAME CS::Vertex_handle new_vertex_cached(CGAL_AOS3_TYPENAME CS::Point p);
  CGAL_AOS3_TYPENAME CS::Halfedge_handle new_halfedge(CGAL_AOS3_TYPENAME CS::Point s, 
						      CGAL_AOS3_TYPENAME CS::Curve ff, 
						      CGAL_AOS3_TYPENAME CS::Point f);
  Combinatorial_cross_section CGAL_AOS3_TARG &cs_;
  const Traits& traits_;
  std::map<Edge, CGAL_AOS3_TYPENAME CS::Halfedge_handle> unmatched_hedges_;
  std::map<CGAL_AOS3_TYPENAME CS::Point, CGAL_AOS3_TYPENAME CS::Vertex_handle> points_;
};
CGAL_AOS3_END_INTERNAL_NAMESPACE


#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Cross_section_initializer_impl.h>
#endif

#endif
