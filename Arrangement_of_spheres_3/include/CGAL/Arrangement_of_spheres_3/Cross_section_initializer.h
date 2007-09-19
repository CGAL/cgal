#ifndef CGAL_ARRANGEMENT_OF_SPHERES_CROSS_SECTION_INIT_H
#define CGAL_ARRANGEMENT_OF_SPHERES_CROSS_SECTION_INIT_H

#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_arrangement.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
struct Cross_section_initializer {
  CGAL_AOS3_TRAITS;
public:
  typedef Combinatorial_cross_section CGAL_AOS3_TARG CS;
  typedef CGAL_AOS3_TYPENAME CS::Point Point;
  typedef CGAL_AOS3_TYPENAME CS::Curve Curve;
  typedef CGAL_AOS3_TYPENAME CS::Halfedge_handle Halfedge_handle;
  typedef CGAL_AOS3_TYPENAME CS::Vertex_handle Vertex_handle;
  typedef CGAL_AOS3_TYPENAME CS::Face_handle Face_handle;
  typedef CGAL_AOS3_TYPENAME Traits::Sphere_3_key Sphere_3_key;
  typedef Cross_section_arrangement CGAL_AOS3_TARG Arr;
  Cross_section_initializer (CS& cs, const Traits& tr);

  void operator()(CGAL_AOS3_TYPENAME Traits::FT z);

  /*struct Edge {
    
    Edge(Point s, 
	 Curve su,
	 Point t): sup_(su), s_(s), t_(t){}
    
    bool operator<(const Edge &o) const {
      if (sup_ < o.sup_) return true;
      else if (o.sup_ < sup_) return false;
      else if (s_ < o.s_) return true;
      else if (o.s_ < s_) return false;
      else return t_ < o.t_;
    }
    
    Curve sup_;
    Point s_, t_;
    };*/
  

  template <class It>
  Halfedge_handle new_face(It b, It e, const Arr &arr) {
    // vt must be a pair, ick 
    CGAL_assertion(std::distance(b,e)>=2);
    It em1=e;
    --em1;
    //CGAL_AOS3_TYPENAME CS::Point lp= arr.point(em1->first);
   
    for (It c=b; c!= e; ++c) {
      /*Halfedge_handle h=*/find_halfedge(*c, arr);
      //CGAL_LOG_WRITE(Log::SOME, cs_.write(h, LOG_STREAM) << ": ");
    }
    CGAL_LOG(Log::SOME, "Face is ");
    for (It c=b; c!= e; ++c) {
      Halfedge_handle h= find_halfedge(*c, arr);
      CGAL_LOG_WRITE(Log::SOME, cs_.write(h, LOG_STREAM) << ": ");
    }
    CGAL_LOG(Log::SOME, std::endl);
      
    CGAL_AOS3_TYPENAME CS::Face_handle f= cs_.hds_.faces_push_back(CGAL_AOS3_TYPENAME CS::HDS::Face());
    It c=b;
    CGAL_AOS3_TYPENAME CS::Halfedge_handle le, fst;
    while (c != e) {
      //CGAL_AOS3_TYPENAME CS::Point cp= arr.point(c->first);
      CGAL_AOS3_TYPENAME Arr::Halfedge_const_handle hh=*c;
      Halfedge_handle h= find_halfedge(hh, arr);
      if (c==b) fst=h;
      h->set_face(f);
      f->set_halfedge(h);
      //lp=cp;
      if (le != CGAL_AOS3_TYPENAME CS::Halfedge_handle()) {
	cs_.connect(le, h);
	//le->set_next(h);
	//h->set_prev(le);
	CGAL_assertion(le->next()->opposite()->vertex() == le->vertex());
      }
      ++c;
      le=h;
    }
    le->set_next(fst);
    fst->set_prev(le);
    return fst;
  }

  ~Cross_section_initializer();
private:
  

  void copy_in(const Arr &arr,
	       std::vector<Halfedge_handle> &ties) ;

  Vertex_handle find_vertex(CGAL_AOS3_TYPENAME Arr::Vertex_const_handle p, const Arr &);


  Halfedge_handle find_halfedge(CGAL_AOS3_TYPENAME Arr::Halfedge_const_handle, const Arr &);


  Combinatorial_cross_section CGAL_AOS3_TARG &cs_;
  const Traits& traits_;
  std::map<CGAL_AOS3_TYPENAME Arr::Halfedge_const_handle,
	   Halfedge_handle,
	   CGAL_AOS3_TYPENAME CS::Handle_compare> unmatched_hedges_;
  std::map<CGAL_AOS3_TYPENAME Arr::Vertex_const_handle,
	   Vertex_handle, CGAL_AOS3_TYPENAME CS::Handle_compare> points_;
};
CGAL_AOS3_END_INTERNAL_NAMESPACE


#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Cross_section_initializer_impl.h>
#endif

#endif
