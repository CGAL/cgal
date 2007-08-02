#ifndef CGAL_ARRANGEMENT_OF_SPHERES_CROSS_SECTION_H
#define CGAL_ARRANGEMENT_OF_SPHERES_CROSS_SECTION_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_face_base.h>
#include <CGAL/Arrangement_of_spheres_3/Rule_direction.h>

#include <CGAL/Arrangement_of_spheres_3/Combinatorial_vertex.h>
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_curve.h>
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <boost/utility.hpp>
#include <map>
#include <set>
#include <boost/array.hpp>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Cross_section_initializer;

CGAL_AOS3_TEMPLATE
class Combinatorial_cross_section: public boost::noncopyable {
public:
  CGAL_AOS3_TRAITS;

  friend class Cross_section_initializer CGAL_AOS3_TARG;

  /*
    Use an hds. For each halfedge mark what sphere, what part and if
    part of the circle, whether inside or outside

    For each edge we also might have an event
  */
public:
  typedef Combinatorial_vertex Point;
  typedef Combinatorial_curve Curve;
  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Rule_direction Rule_direction;
 
  typedef CGAL_AOS3_TYPENAME Traits::Event_key Event_key;

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
	Event_key event() const {
	  return ev_;
	}
	void set_event(Event_key ev) {
	  ev_=ev;
	}

	Event_key ev_;
	Curve pt_;
      };
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
      // maybe want pointers to chains, maybe not
      typedef CGAL::HalfedgeDS_face_base< Refs> Face;
    };
  };
  

    

  typedef CGAL::HalfedgeDS_default<int, Slice_halfedgeDS_items_2> HDS;

  typedef CGAL_AOS3_TYPENAME HDS::Halfedge_handle Halfedge_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Halfedge_const_handle Halfedge_const_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Vertex_handle Vertex_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Vertex_const_handle Vertex_const_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Face_handle Face_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Face_const_handle Face_const_handle;

  Combinatorial_cross_section();

  HDS& hds();

  struct HE_key: public HDS::Halfedge_const_handle {
    typedef CGAL_AOS3_TYPENAME HDS::Halfedge_const_handle P;
    HE_key(P h): P(h){}
    bool operator<(P o) const {
      return &P::operator*() < &o.operator*();
    }
  };



  CGAL_CONST_ITERATOR(Halfedge, halfedge, CGAL_AOS3_TYPENAME HDS::Halfedge_const_iterator,
		      return hds_.halfedges_begin(),
		      return hds_.halfedges_end());
  
  
  CGAL_CONST_ITERATOR(Vertex, vertice, CGAL_AOS3_TYPENAME HDS::Vertex_const_iterator,
		      return hds_.vertices_begin(),
		      return hds_.vertices_end());
  
  CGAL_CONST_ITERATOR(Face, face, CGAL_AOS3_TYPENAME HDS::Face_const_iterator,
		      return hds_.faces_begin(),
		      return hds_.faces_end());
  
  CGAL_SIZE(vertices, return hds_.size_of_vertices();)

  CGAL_SIZE(faces, return hds_.size_of_faces();)
  CGAL_SIZE(edges, return hds_.size_of_halfedges();)

  Halfedge_handle find_halfedge(Vertex_handle v, Face_handle f);

  bool has_vertex(Face_const_handle fh, Vertex_const_handle vh) const;

  std::ostream &write(Halfedge_const_handle h, std::ostream &out) const;

  std::ostream &write(Face_const_handle h, std::ostream &out) const;

  std::ostream &write(std::ostream &out) const;

  bool is_in_slice(Vertex_const_handle v) const;
  bool is_in_slice(Halfedge_const_handle h) const;
  bool is_in_slice(Face_const_handle h) const;

  unsigned int degree(Vertex_const_handle v) const;

  bool is_redundant(Vertex_const_handle v) const;

  void clear();

  CGAL_GETNR(Face_handle, infinite_face, return inf_);
  //Face_handle infinite_face(){return inf_;}

  CGAL_ITERATOR(Halfedge, halfedge, CGAL_AOS3_TYPENAME HDS::Halfedge_iterator,
		return hds_.halfedges_begin(),
		return hds_.halfedges_end());

  CGAL_ITERATOR(Face, face,CGAL_AOS3_TYPENAME HDS::Face_iterator,
		return hds_.faces_begin(),
		return hds_.faces_end());

  
  // the first is the target of h, second the source
  std::pair<Halfedge_handle, Halfedge_handle> remove_edge(Halfedge_handle h);

  void set_curve(Halfedge_handle h, Curve c);
   
  /*void set_event(Halfedge_handle h, Event_key k) const {
    CGAL_assertion(k != Event_key());
    CGAL_assertion(h->event()== h->opposite()->event());
    CGAL_assertion(h->event() == Event_key());
    h->set_event(k);
    h->opposite()->set_event(k);
  }

  void unset_event(Halfedge_handle h) const {
    CGAL_assertion(h->event()== h->opposite()->event());
    h->set_event(Event_key());
    h->opposite()->set_event(Event_key());
    }*/

  Face_handle merge_faces(Halfedge_handle e);

  Face_handle merge_faces(Vertex_handle v);
  
  Halfedge_handle split_face(Curve c, Halfedge_handle source,
			     Halfedge_handle target);
  
  void insert_target(Curve::Key k,
		     Halfedge_handle cps[]);

  Face_handle remove_target(Halfedge_handle ts[],
			    Halfedge_handle vert[]);
  /*// move old_edge to have the vertices be in new_source and new_target
    // return the new edge
    Halfedge_handle move_edge(Halfedge_handle old_edge,
    Halfedge_handle new_source,
    Halfedge_handle new_target,
    Halfedge_handle &blank_old_source,
    Halfedge_handle &blank_old_target);*/
	
  // remove the rule from the endpoint
  // attach it to the new one
  void move_edge_target(Halfedge_handle edge,
			Halfedge_handle new_target); 

  //! stitch the the HDS together by attaching the vertices pointed to by ib:ie to those in fb,fe
  /*!
    New edges are stuck in edges
   */
  template <class It, class Cit, class Oit>
  void stitch_in(It ib, It ie, It fb, Cit cb, Oit edges) {
    Face_handle f= (*fb)->begin();
    for (It ic=ib; ic != ie; ++ic) {
      CGAL_assertion((*ic)->face() == void_);
    }
    Halfedge_handle c= *ib;
    do {
      c->set_face(f);
      c=c->next();
    } while (c != *ib);
    
    Halfedge_handle hn= new_halfedge(*cb);
    *edges= hn;
    ++edges;
    hn->set_face(f);
    hn->opposite()->set_face(f);
    connect( hn, (*fb)->next());
    connect(*fb, hn->opposite());
    connect(*ib, hn);
    connect(hn->opposite(), (*ib)->next());

    ++ib; ++fb;
    CGAL::HalfedgeDS_decorator<HDS> hdsd(hds_);
    for (; ib != ie; ++ib) {
      Halfedge_handle h= hdsd.split_face(*ib, *fb);
      h->set_curve(*cb);
      h->opposite()->set_curve(cb->other_side());
      *edges=h;
      ++edges;
    }
  }

  void audit() const ;	 


  typedef std::pair<Point, Curve> NFP;

  void reserve(int nv, int ne, int nf);

 
  void set_has_boundary(bool tf);

  bool has_boundary() const;

  Halfedge_handle remove_vertex(Halfedge_handle v);

  void exchange_sphere_extremums(Curve::Key k, Curve::Key l);

  Halfedge_handle halfedge(Vertex_handle v, Face_handle f) const {
    Halfedge_handle h= v->halfedge();
    do {
      if (h->face() == f) return h;
      h= h->next()->opposite();
    } while(h != v->halfedge());
    CGAL_assertion(0);
    return Halfedge_handle();
  }

  Halfedge_handle next_edge_on_curve(Halfedge_handle) const;

  Halfedge_handle cross_edge(Halfedge_handle) const;

  // a halfedge on the curve (an inside one)
  Halfedge_handle a_halfedge(Curve::Key k) const;

  // a halfedge on the rule pointing to the extremal vertex
  Halfedge_handle rule_halfedge(Curve::Key k, Rule_direction i) const;

  // a halfedge on the circle pointing to the extremal vertex (inside)
  Halfedge_handle extremum_halfedge(Curve::Key k, Rule_direction i) const;
  
  // insert the vertex so that h->opposite points to it
  Halfedge_handle insert_vertex( Point p, Halfedge_handle h);


  std::pair<Halfedge_handle,Halfedge_handle> 
  intersect(Halfedge_handle ha, Halfedge_handle hb);

  std::pair<Halfedge_handle,Halfedge_handle>  unintersect(Face_handle f);

 
  
  //friend class Cross_section_initializer CGAL_AOS3_TARG;
 

  void relabel_rule(Halfedge_handle h, Curve nl);

 

  typedef boost::array<Halfedge_handle, 4> Halfedge_quad;

  const Halfedge_quad& halfedges(CGAL_AOS3_TYPENAME Curve::Key k) const {
    return halfedges_[k.input_index()];
  }

   Halfedge_quad& halfedges(CGAL_AOS3_TYPENAME Curve::Key k)  {
    return halfedges_[k.input_index()];
  }
  void set_number_of_spheres(unsigned int i) {
    CGAL_assertion(halfedges_.size() <=i);
    halfedges_.resize(i);
  }


private:
  void connect(Halfedge_handle a, Halfedge_handle b);

  void new_circle(Curve::Key k, Face_handle f, Halfedge_handle c[]);
  
  /*typedef boost::tuple<Halfedge_handle, Halfedge_handle,
    Halfedge_handle, Halfedge_handle> Halfedge_handle_quadruple;*/

  void new_target(Curve::Key k, Halfedge_handle tar[]);
  
  void set_extremum_halfedge(Halfedge_handle h);


  /*Halfedge_handle split_face(Halfedge_handle o, Halfedge_handle d,
    Curve c);*/

  void relabel_target(Halfedge_handle v[], Curve::Key k);

  
  //void exchange_vertices(Halfedge_handle h, Halfedge_handle p);

  void audit_vertex(Vertex_const_handle v) const ;
  


  //void audit_halfedge(Halfedge_const_handle v) const ;





  void delete_edge(Halfedge_handle h);

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
  Halfedge_handle insert_vertex_in_edge_unsafe(Halfedge_handle h,
					       Vertex_handle  vh);

  // insert a vertex for p in both edges
  // merge the two vertices
  std::pair<Halfedge_handle, Halfedge_handle> pinch_bl(Halfedge_handle a,
						       Halfedge_handle b, Point p);

  // a and b go to the same vertex. Make them go to different vertices
  // return the new vertex on b
  Vertex_handle unpinch_bl(Halfedge_handle a, Halfedge_handle b);

  void merge_vertices_bl(Halfedge_handle a, Halfedge_handle b);

  Halfedge_handle new_halfedge(Curve c);

  //typedef std::pair<Point,Point>  ED;

  Vertex_handle new_vertex(Point p);


 
protected:


  mutable HDS hds_;
  Face_handle inf_;
  Face_handle void_;
  mutable std::vector<Curve> errors_;
  std::vector<Vertex_handle> targets_;
  std::vector<Halfedge_quad> halfedges_;

  // for new_face when constructing things
  //std::map<Edge, Halfedge_handle> unmatched_hedges_;
  //std::map<Point, Vertex_handle> points_;
    
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section_impl.h>
#endif

#endif
