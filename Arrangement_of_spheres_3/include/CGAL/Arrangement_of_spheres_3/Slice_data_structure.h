#ifndef SLICE_DATA_STRUCTURE_H
#define SLICE_DATA_STRUCTURE_H

#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_face_base.h>

#include <CGAL/Arrangement_of_spheres_3/Simulator.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_vertex.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>
#include <map>
#include <set>
#include <boost/array.hpp>



class Slice_data_structure {
public:

  /*
    Use an hds. For each halfedge mark what sphere, what part and if
    part of the circle, whether inside or outside

    For each edge we also might have an event
  */

  typedef Combinatorial_vertex Point;
  typedef Combinatorial_curve Curve;


  struct Slice_halfedgeDS_items_2 {
   

    template < class Refs, class Traits>
    struct Vertex_wrapper {
     
      struct Vertex: public CGAL::HalfedgeDS_vertex_base< Refs, CGAL::Tag_true, Point>{
      };
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper {
      struct Halfedge: public CGAL::HalfedgeDS_halfedge_base< Refs> {
	Halfedge(){}
	Curve curve() const {
	  return pt_;
	}
	Curve& curve() {
	  return pt_;
	}
	void set_curve(Curve pt) {
	  pt_=pt;
	}
	typename Simulator::Event_key event() const {
	  return ev_;
	}
	void set_event(Simulator::Event_key ev) {
	  ev_=ev;
	}

	Simulator::Event_key ev_;
	Curve pt_;
      };
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
      // maybe want pointers to chains, maybe not
      typedef CGAL::HalfedgeDS_face_base< Refs> Face;
    };
  };
  

  /* struct Next {
    Halfedge_handle operator()(Halfedge_handle h) const {
      return h->next();
    }
  };
  struct Prev {
    Halfedge_handle operator()(Halfedge_handle h) const {
      return h->prev();
    }
    };*/

  typedef CGAL::HalfedgeDS_default<int, Slice_halfedgeDS_items_2> HDS;

  typedef HDS::Halfedge_handle Halfedge_handle;
  typedef HDS::Halfedge_const_handle Halfedge_const_handle;
  typedef HDS::Vertex_handle Vertex_handle;
  typedef HDS::Vertex_const_handle Vertex_const_handle;
  typedef HDS::Face_handle Face_handle;
  typedef HDS::Face_const_handle Face_const_handle;

  Slice_data_structure(int num);

  HDS& hds();

  struct HE_key: public HDS::Halfedge_const_handle {
    typedef HDS::Halfedge_const_handle P;
    HE_key(P h): P(h){}
    bool operator<(P o) const {
      return &P::operator*() < &o.operator*();
    }
  };

 
  void debug_add_sphere(){
    halfedges_.resize(halfedges_.size()+1);
  }

  typedef HDS::Halfedge_const_iterator Halfedge_const_iterator;
  Halfedge_const_iterator halfedges_begin() const {
    return hds_.halfedges_begin();
  }
  Halfedge_const_iterator halfedges_end() const {
    return hds_.halfedges_end();
  }


  typedef HDS::Halfedge_iterator Halfedge_iterator;
  Halfedge_iterator halfedges_begin() {
    return hds_.halfedges_begin();
  }
  Halfedge_iterator halfedges_end() {
    return hds_.halfedges_end();
  }

  typedef HDS::Vertex_const_iterator Vertex_const_iterator;
  Vertex_const_iterator  vertices_begin() const {
    return hds_.vertices_begin();
  }
  Vertex_const_iterator vertices_end() const {
    return hds_.vertices_end();
  }

  typedef HDS::Face_const_iterator Face_const_iterator;
  Face_const_iterator  faces_begin() const {
    return hds_.faces_begin();
  }
  Face_const_iterator faces_end() const {
    return hds_.faces_end();
  }

  typedef HDS::Face_iterator Face_iterator;
  Face_iterator  faces_begin() {
    return hds_.faces_begin();
  }
  Face_iterator faces_end() {
    return hds_.faces_end();
  }

  Halfedge_handle find_halfedge(Vertex_handle v, Face_handle f);

  // the first is the target of h, second the source
  std::pair<Halfedge_handle, Halfedge_handle> remove_rule(Halfedge_handle h);

  bool has_vertex(Face_const_handle fh, Vertex_const_handle vh) const;

  void connect(Halfedge_handle a, Halfedge_handle b);

  Simulator::Event_key event(Halfedge_handle h) const {
    return h->event();
  }

  void set_event(Halfedge_handle h, Simulator::Event_key k) const {
    CGAL_assertion(k != Simulator::Event_key());
    CGAL_assertion(h->event()== h->opposite()->event());
    CGAL_assertion(h->event() == Simulator::Event_key());
    h->set_event(k);
    h->opposite()->set_event(k);
  }

  void unset_event(Halfedge_handle h) const {
    CGAL_assertion(h->event()== h->opposite()->event());
    h->set_event(Simulator::Event_key());
    h->opposite()->set_event(Simulator::Event_key());
  }

  void new_circle(Curve::Key k, Face_handle f, Halfedge_handle c[]);
  
  /*typedef boost::tuple<Halfedge_handle, Halfedge_handle,
    Halfedge_handle, Halfedge_handle> Halfedge_handle_quadruple;*/

  void new_target(Curve::Key k, Halfedge_handle tar[]);
  
  void insert_target(Curve::Key k,
		     Halfedge_handle cps[]);

  Face_handle remove_target(Halfedge_handle ts[],
			    Halfedge_handle vert[]);

  Halfedge_handle split_face(Halfedge_handle o, Halfedge_handle d,
			     Curve c);

  void relabel_target(Halfedge_handle v[], Curve::Key k);

  
  void exchange_vertices(Halfedge_handle h, Halfedge_handle p);

  Face_handle intersect(Halfedge_handle ha, Halfedge_handle hb);

  std::pair<Halfedge_handle,Halfedge_handle>  unintersect(Face_handle f);

  // remove the rule from where it is connected
  // attach it to point a a new vertex in ne, starting from t->vertex()
  // the old target is removed, the old source is not. 
  // OK, this is really screwy
  // the new rule is returned
  void move_rule(Halfedge_handle r,
			    Halfedge_handle ne, 
			    Halfedge_handle t); 


  void audit() const ;

  void audit_vertex(Vertex_const_handle v) const ;
  
  Halfedge_handle next_edge_on_curve(Halfedge_handle) const;
  Halfedge_const_handle next_edge_on_curve(Halfedge_const_handle) const;
  Halfedge_handle cross_edge(Halfedge_handle) const;

  //void audit_halfedge(Halfedge_const_handle v) const ;

  void set_is_building(bool tf);
  
  std::ostream &write(Halfedge_const_handle h, std::ostream &out) const;

  std::ostream &write_face(Halfedge_const_handle h, std::ostream &out) const;

  bool is_in_slice(Vertex_const_handle v) const;
  bool is_in_slice(Halfedge_const_handle h) const;
  bool is_in_slice(Face_const_handle h) const;

  typedef std::pair<Point, Curve> NFP;

  void reserve(int nv, int ne, int nf);

  void initialize(int num) ;

  unsigned int degree(Vertex_const_handle v) const;

  bool is_redundant(Vertex_const_handle v) const;

  Halfedge_handle remove_redundant_vertex(Halfedge_handle v);

  void clear();

  // a halfedge on the curve (an inside one)
  Halfedge_handle a_halfedge(Curve::Key k) const;

  // a halfedge on the rule pointing to the extremal vertex
  Halfedge_handle rule_halfedge(Curve::Key k, int i) const;


  template <class Out>
  void find_halfedges(Curve c, Out o) const {
    for (Halfedge_const_iterator it= halfedges_begin();
	 it != halfedges_end(); ++it) {
      if (it->curve() == c) {
	*o=it;
	++o;
      }
    }
  }


  // insert the vertex so that h->opposite points to it
  Vertex_handle insert_vertex_in_edge(Halfedge_handle h, Point p);

  Halfedge_handle new_halfedge(Curve c);

  //typedef std::pair<Point,Point>  ED;

  Vertex_handle new_vertex(Point p);

  Vertex_handle new_vertex_cached(Point p);

  Halfedge_handle new_halfedge(Point s, Curve ff, Point f);

  template <class It>
  void new_face(It b, It e) {
    // vt must be a pair, ick 
    CGAL_assertion(std::distance(b,e)>=2);
    It em1=e;
    --em1;
    Point lp= em1->second;
    if (points_.find(lp)==points_.end()) new_vertex_cached(lp);
   
    Face_handle f= hds_.faces_push_back(HDS::Face());
    It c=b;
    Halfedge_handle le, fst;
    while (c != e) {
      Point cp= c->second;
      Halfedge_handle h= new_halfedge(lp, c->first, cp);
      if (c==b) fst=h;
      h->set_face(f);
      f->set_halfedge(h);
      lp=cp;
      if (le != Halfedge_handle()) {
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


  struct Edge {
    
    Edge(Point s, Curve su, Point t): sup_(su), s_(s), t_(t){}
    
    bool operator<(const Edge &o) const {
      if (sup_ < o.sup_) return true;
      else if (o.sup_ < sup_) return false;
      else if (s_ < o.s_) return true;
      else if (o.s_ < s_) return false;
      else return t_ < o.t_;
    }
    
    Curve sup_;
    Point s_, t_;
  };
  

  mutable HDS hds_;
  Face_handle inf_;
  mutable std::vector<Curve> errors_;
  std::vector<Vertex_handle> targets_;
  typedef boost::array<Halfedge_handle, 4> Halfedge_quad;
  std::vector<Halfedge_quad> halfedges_;

  // for new_face when constructing things
  std::map<Edge, Halfedge_handle> unmatched_hedges_;
  std::map<Point, Vertex_handle> points_;
};

#endif
