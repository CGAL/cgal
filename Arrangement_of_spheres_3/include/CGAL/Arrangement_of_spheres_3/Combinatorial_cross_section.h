#ifndef CGAL_ARRANGEMENT_OF_SPHERES_CROSS_SECTION_H
#define CGAL_ARRANGEMENT_OF_SPHERES_CROSS_SECTION_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/HalfedgeDS_default.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_decorator.h>
#include <CGAL/HalfedgeDS_face_base.h>
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Arrangement_of_spheres_3/Event_visitor.h>
#include <CGAL/Arrangement_of_spheres_3/Cross_section_halfedgeDS_items_2.h>
#include <boost/utility.hpp>
#include <map>
#include <set>
//#include <boost/array.hpp>
#include <boost/foreach.hpp>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Cross_section_initializer;

CGAL_AOS3_TEMPLATE
class Combinatorial_cross_section: public boost::noncopyable {
public:
  CGAL_AOS3_TRAITS;
  typedef Combinatorial_cross_section CGAL_AOS3_TARG This;
public:
  friend class Cross_section_initializer CGAL_AOS3_TARG;

  /*
    Use an hds. For each halfedge mark what sphere, what part and if
    part of the circle, whether inside or outside

    For each edge we also might have an event
  */

  struct Handle_compare{
    template <class H>
    bool operator()(H a, H b) const {
      return &*a < &*b;
    }
  };

  struct Handle_equal{
    template <class H>
    bool operator()(H a, H b) const {
      return &*a == &*b;
    }
  };

  

  typedef Combinatorial_vertex Point;
  typedef Combinatorial_curve Curve;
  typedef CGAL_AOS3_TYPENAME CGAL_AOS3_INTERNAL_NS::Rule_direction Rule_direction;
 
  typedef CGAL_AOS3_TYPENAME Traits::Event_key Event_key;

  

    

  typedef CGAL::HalfedgeDS_default<int, Cross_section_halfedgeDS_items_2 CGAL_AOS3_TARG > HDS;

  typedef CGAL_AOS3_TYPENAME HDS::Halfedge_handle Halfedge_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Halfedge_const_handle Halfedge_const_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Vertex_handle Vertex_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Vertex_const_handle Vertex_const_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Face_handle Face_handle;
  typedef CGAL_AOS3_TYPENAME HDS::Face_const_handle Face_const_handle;

  typedef Event_visitor CGAL_AOS3_TARG Visitor;



  Combinatorial_cross_section();

  HDS& hds();

  struct HE_key: public HDS::Halfedge_const_handle {
    typedef CGAL_AOS3_TYPENAME HDS::Halfedge_const_handle P;
    HE_key(P h): P(h){}
    bool operator<(P o) const {
      return &P::operator*() < &o.operator*();
    }
  };


  Visitor& visitor() {
    return v_;
  }


  const Visitor& visitor()const {
    return v_;
  }

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
  std::ostream &write(Vertex_const_handle h, std::ostream &out) const;

  std::ostream &write(Face_const_handle h, std::ostream &out) const;

  std::ostream &write(std::ostream &out) const;

  bool is_in_slice(Vertex_const_handle v) const;
  bool is_in_slice(Halfedge_const_handle h) const;
  bool is_in_slice(Face_const_handle h) const;

  unsigned int degree(Vertex_const_handle v) const;

  bool is_redundant(Vertex_const_handle v) const;
  bool is_redundant(Halfedge_const_handle v) const;

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

  void init_halfedge(Halfedge_handle h, Curve c);
   
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

  Face_handle join_face(Halfedge_handle e, bool cleanup);

  //Face_handle join_face(Vertex_handle v);
  
  Halfedge_handle split_face(Curve c, Halfedge_handle source,
			     Halfedge_handle target);
  
  /*void insert_target(Curve::Key k,
    Halfedge_handle cps[]);

    Face_handle remove_target(Halfedge_handle ts[],
    Halfedge_handle vert[]);*/
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

 
  template <class Hit, class Cit>
  Vertex_handle star_face(const Hit hb, Hit he, const Cit cb, Cit ce) {
    CGAL_LOG(Log::SOME, "Star face called for face");
    CGAL_LOG_WRITE(Log::SOME, write((*hb)->face(), LOG_STREAM) << std::endl);
    CGAL_LOG(Log::SOME, "Handles are ");
    {
      Cit cc=cb;
      for (Hit hc=hb; hc != he; ++hc, ++cc) {
	CGAL_LOG_WRITE(Log::SOME, write(*hc, LOG_STREAM) 
		       << ": " << *cc << std::endl);
      }
    }
    CGAL_assertion(std::distance(hb,he) == std::distance(cb, ce));
    Vertex_handle vh= new_vertex(Point());
    Face_handle f= (*hb)->face();
    for (Hit hc=hb; hc != he; ++hc) {
      CGAL_assertion((*hc)->face() == f);
    }
    Hit hc=hb;
    Cit cc=cb;
    Halfedge_handle hf= new_halfedge(cc->opposite());
    hf->set_face(f);
    hf->opposite()->set_face(f);
    hf->set_vertex(vh);
    vh->set_halfedge(hf);
    hf->opposite()->set_vertex((*hc)->vertex());
    connect(hf->opposite(), (*hc)->next());
    connect(*hc, hf);
    connect(hf, hf->opposite());
    ++hc; ++cc;

    CGAL::HalfedgeDS_decorator<HDS> hdsd(hds_);
    //CGAL::HalfedgeDS_const_decorator<HDS> hdscd(hds_);
    CGAL_LOG_WRITE(Log::SOME, write(hf->face(), LOG_STREAM) << " is the new face ");
    //hdscd.is_valid(true, 3);
    for (; hc != he; ++hc, ++cc) {
      CGAL_assertion(hf->face() == (*hc) ->face());
      CGAL_LOG(Log::SOME, "Connecting ");
      CGAL_LOG_WRITE(Log::SOME, write(hf, LOG_STREAM) << " to ");
      CGAL_LOG_WRITE(Log::SOME, write(*hc, LOG_STREAM) << std::endl);
      Halfedge_handle hn=hdsd.split_face(hf, *hc);
      set_curve(hn, *cc);
    }
    {
      CGAL_LOG(Log::SOME, "Result is\n");
      Halfedge_handle h= vh->halfedge();
      do {
	CGAL_LOG_WRITE(Log::SOME, write(h, LOG_STREAM) << std::endl);
	h= h->next()->opposite();
      } while (h != vh->halfedge());
    }
    return vh;
  }


  template <class Iit, class Oit>  
  void expand_vertex(Iit ib, Iit ie,
		     Oit ob, Oit oe) {
    CGAL_assertion(std::distance(ib,ie) == std::distance(ob, oe));
    CGAL_LOG(Log::SOME, "Explanding vertex with\n");
    {
      Oit oc=ob;
      for (Iit ic=ib; ic != ie; ++ic, ++oc) {
	CGAL_LOG_WRITE(Log::SOME, write(*ic, LOG_STREAM) << ": ");
	CGAL_LOG_WRITE(Log::SOME, write(*oc, LOG_STREAM) << std::endl);
      }
    }
    Vertex_handle vh= (*ob)->vertex();
    Face_handle f=(*ib)->face();
    CGAL_assertion(f != inf_);
    {
      Oit oc=ob;
      for (Iit ic=ib; ic != ie; ++ic, ++oc) {
	CGAL_LOG_WRITE(Log::SOME, write((*ic)->opposite(), LOG_STREAM) << ": ");
	CGAL_LOG_WRITE(Log::SOME, write(*oc, LOG_STREAM) << std::endl);
	CGAL_assertion((*ic)->opposite()->curve() == (*oc)->curve());
	CGAL_assertion(degree((*ic)->vertex()) ==1);
	CGAL_assertion((*oc)->vertex() ==vh);
	CGAL_assertion((*ic)->face() ==f);
      }
    }
    

    CGAL_LOG(Log::LOTS,  "Auditing const decorator..." << std::flush);
    CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
    if (!chds.is_valid(false, num_components_==1? 3: 2)) {
      chds.is_valid(true, num_components_==1? 3: 2);
      std::cerr << "Not valid." << std::endl;
      CGAL_error();
    }

    {
      Oit oc=ob;
      for (Iit ic=ib; ic != ie; ++ic, ++oc) {
	Halfedge_handle hi= (*ic)->opposite();
	Halfedge_handle ho= *oc;
	CGAL_assertion(hi->curve() == ho->curve());
	connect(ho, hi->next());
	connect(hi->opposite()->prev(), ho->opposite());
	ho->set_vertex(hi->vertex());
	ho->vertex()->set_halfedge(ho);
      }
    }
    {
      CGAL::HalfedgeDS_decorator<HDS> hdsd(hds_);
      Oit oc=ob;
      for (Iit ic=ib; ic != ie; ++ic, ++oc) {
	hdsd.set_face_in_face_loop(*oc, (*oc)->face());
	hds_.vertices_erase((*ic)->vertex());
	hds_.edges_erase(*ic);
      }
    }
    hds_.vertices_erase(vh);
    hds_.faces_erase(f);
    --num_components_;
    //std::set<Halfedge_handle, Handle_compare> visited;
    for (Oit oc=ob; oc != oe; ++oc) {
      v_.on_new_edge(*oc);
    }
  }

  //! stitch the the HDS together by attaching the vertices pointed to by ib:ie to those in fb,fe
  /*!
    New edges are stuck in edges
  */
  template <class It, class Cit>
  void stitch_in(It ib, It ie, It fb, Cit cb) {
    CGAL_precondition(num_components_ > 1);
    CGAL_LOG(Log::LOTS, "Stitching: " << std::endl);
    Cit cc=cb;
    for (It ic= ib, fc= fb; ic != ie; ++ic, ++fc, ++cc) {
      CGAL_LOG_WRITE(Log::LOTS, write(*ic, LOG_STREAM) << " to ");
      CGAL_LOG_WRITE(Log::LOTS, write(*fc, LOG_STREAM) << " with " << *cc << std::endl);
    }

    for (It ic= ib, fc= fb; ic != ie; ++ic, ++fc, ++cc) {
      CGAL_assertion(is_valid(*ic));
      CGAL_assertion(is_valid(*fc));
    }

    --num_components_;
    Face_handle f= (*fb)->face();
    for (It ic=ib, fc= fb; ic != ie; ++ic, ++fc) {
      CGAL_assertion((*ic)->face() == inf_);
      CGAL_assertion((*fc)->face() == (*fb)->face());
    }
   
    Halfedge_handle c= *ib;
    do {
      CGAL_assertion(c->face() == inf_);
      c->set_face(f);
      c=c->next();
    } while (c != *ib);
    
    Halfedge_handle hn= new_halfedge(*cb);
    //*edges= hn;
    //++edges;
    hn->set_face(f);
    hn->opposite()->set_face(f);
    hn->set_vertex((*fb)->vertex());
    hn->opposite()->set_vertex((*ib)->vertex());
    connect( hn, (*fb)->next());
    connect(*fb, hn->opposite());
    connect(hn->opposite(), (*ib)->next());
    connect(*ib, hn);

    std::vector<Halfedge_handle> dirty;

    dirty.push_back(hn);
    v_.on_delete_edge(hn->next());
    v_.on_delete_edge(hn->prev());
    v_.on_delete_edge(hn->opposite()->next());
    v_.on_delete_edge(hn->opposite()->prev());
    dirty.push_back(hn->next());
    dirty.push_back(hn->prev());
    dirty.push_back(hn->opposite()->next());
    dirty.push_back(hn->opposite()->prev());
    		  

    ++ib; ++fb; ++cb;
    CGAL::HalfedgeDS_decorator<HDS> hdsd(hds_);
   
    for (; ib != ie; ++ib, ++fb, ++cb) {
      //split_face(*cb, *ib, *fb);
      v_.on_delete_edge(*ib);
      v_.on_delete_edge((*ib)->next());
      dirty.push_back(*ib);
      dirty.push_back((*ib)->next());
      v_.on_delete_edge(*fb);
      v_.on_delete_edge((*fb)->next());
      dirty.push_back(*fb);
      dirty.push_back((*fb)->next());
      Halfedge_handle h= hdsd.split_face(*ib, *fb);
      h->set_curve(*cb);
      h->opposite()->set_curve(cb->opposite());
      
      dirty.push_back(h);
      
      
      // *edges=h;
      //++edges;
    }
    handle_new_edges(dirty);

    audit();
  }

  // edges must be CCW around hole.

  template <class It>
  Face_handle snip_out(It b, It e, bool clean) {
    // want to make this clean up after itself and properly visit changed edges after everything is done. 

#ifndef NDEBUG
    // edges must be ccw and oriented out
    for (It c= b; c !=e ; ++c){
      It cp1= c;
      ++cp1;
      if (cp1==e) cp1=b;
      Halfedge_handle h= *c;
      h=h->next();
      do {
	for (It ic= b; ic !=e ; ++ic){
	  if (ic == cp1) continue;
	  CGAL_assertion(h != *ic);
	  CGAL_assertion(h->opposite()!= *ic);
	}
	h=h->next();
      } while (h->opposite() != *cp1);
    }
    for (It c= b; c !=e ; ++c){
      CGAL_assertion(is_valid(*c));
    }
#endif
    
    CGAL_LOG(Log::LOTS, "Snipping out ");
    for (It c= b; c != e; ++c) {
      CGAL_LOG_WRITE(Log::LOTS, write(*c, LOG_STREAM) << ", ");
    }
    CGAL_LOG(Log::LOTS, std::endl);
    ++num_components_;
    It c=b;
    ++c; ++c;
    CGAL::HalfedgeDS_decorator<HDS> hdsd(hds_);
    //std::vector<Halfedge_handle> dirty;
    for (; c!= e; ++c) {
      v_.on_delete_edge(*c);
      v_.on_merge_faces(*c);
      Halfedge_handle n= (*c)->next();
      //v_.on_delete_edge(n);
      //v_.on_delete_edge(n->prev());

      //Halfedge_handle p= (*c)->prev();
      hdsd.join_face(*c);
      //dirty.push_back(n);
      //dirty.push_back(n->prev());
     
      //v_.on_changed_edge(p);
      //v_.on_changed_edge(p->next());
    }
    c=b;
    Halfedge_handle h1= *c;
    ++c;
    Halfedge_handle h2= *c;

    Halfedge_handle h1n= h1->next();
    Halfedge_handle h2n= h2->next();
    if (h1n == h2->opposite()) h1n= Halfedge_handle();
    else {
      //v_.on_delete_edge(h1n);
      //dirty.push_back(h1n);
    }
    if (h2n == h1->opposite()) h2n= Halfedge_handle();
    else {
      //v_.on_delete_edge(h2n);
      //dirty.push_back(h2n);
    }


    v_.on_delete_edge(h1);
    v_.on_delete_edge(h2);
    v_.on_merge_faces(h1);

    Halfedge_handle inside, outside;
    Face_handle f= h1->face();
    Face_handle of= h1->opposite()->face();


    // h1,h2 point out
    if (h1->prev() != h2->opposite()) {
      inside = h1->prev();
    } else if (h2->prev() != h1->opposite()) {
      inside= h2->prev();
    } else {
      CGAL_LOG(Log::LOTS, "Single vertex remaining" << std::endl);
      h1->opposite()->vertex()->set_halfedge(Halfedge_handle());
    }
    if (h1->next() != h2->opposite()) {
      outside= h1->next();
    } else if (h2->next() != h1->opposite()) {
      outside= h2->next();
    }

    h1->opposite()->prev()->vertex()->set_halfedge(h1->opposite()->prev());
    h2->opposite()->prev()->vertex()->set_halfedge(h2->opposite()->prev());

    if (inside != Halfedge_handle()) {
      h1->prev()->vertex()->set_halfedge(h1->prev());
      h2->prev()->vertex()->set_halfedge(h2->prev());
    }

    connect(h1->opposite()->prev(), h1->next());
    connect(h2->opposite()->prev(), h2->next());
 
    connect(h1->prev(), h1->opposite()->next());
    connect(h2->prev(), h2->opposite()->next());
   
    CGAL_assertion(outside != Halfedge_handle());
    Halfedge_handle hc= outside;
    do {
      hc->vertex()->set_halfedge(hc);
      hc->set_face(f);
      hc=hc->next();
    } while (hc != h1->next());
    f->set_halfedge(outside);


    if (inside != Halfedge_handle()) {
      Halfedge_handle hc= inside;
      do {
	hc->set_face(inf_);
	hc->vertex()->set_halfedge(hc);
	hc=hc->next();
      } while (hc != inside);
    } else {
      h1->opposite()->vertex()->set_halfedge(Halfedge_handle());
    }
    hds_.faces_erase(of);
    delete_edge(h1);
    delete_edge(h2);
    
    Halfedge_handle ho=h1n;
    if (ho== Halfedge_handle()) ho=h2n;

    if (clean) {
      Halfedge_handle c= outside; 
      do {
	if (is_redundant(c->vertex())) {
	  c=remove_vertex(c->vertex());
	} else {
	  c=c->next();
	}
      } while (c != outside);
    }
    {
      Halfedge_handle c= f->halfedge(); 
      do {
	v_.on_change_edge(c);
	c=c->next();
      } while (c != outside);
    }
    /* // not efficient and not isolated
    if (h1n != Halfedge_handle()) {
      v_.on_change_edge(h1n);
      v_.on_change_edge(h1n->prev());
    }
    if (h2n != Halfedge_handle()) {
      v_.on_change_edge(h2n);
      v_.on_change_edge(h2n->prev());
    }
    */

    write(std::cout);
    /* std::cout << "Auditing const decorator..." << std::flush;
       CGAL::HalfedgeDS_const_decorator<HDS> chds(hds_);
       if (!chds.is_valid(true, num_components_==1? 3: 2)) {
       CGAL_error();
       std::cerr << "Not valid." << std::endl;
       }
       std::cout << "done." << std::endl;*/

    return f;
  }

  void visit_component(Halfedge_handle h) {
    std::set<Halfedge_handle, Handle_compare> visited;
    std::vector<Halfedge_handle> front;
    front.push_back(h); front.push_back(h->opposite());
    do {
      Halfedge_handle h= front.back();
      front.pop_back();
      Halfedge_handle c=h->opposite()->prev();
      while (c != h) {
	if (degree(h->opposite()->vertex()) != 1) {
	  if (visited.find(h->opposite()) != visited.end()) {
	    CGAL_assertion(visited.find(h) == visited.end());
	    v_.on_new_edge(h);
	    visited.insert(h);
	  }
	}
	c= c->opposite()->prev();
      };
    } while (!front.empty());
  }

  Vertex_handle degen_join_vertex(Halfedge_handle h) {
    CGAL::HalfedgeDS_decorator<HDS> hdsd(hds_);
    v_.on_delete_edge(h);
    Vertex_handle vh= h->vertex();
    hdsd.join_vertex(h);
    audit_vertex(vh, false);
    return vh;
  }

  void handle_new_edges(std::vector<Halfedge_handle> &edges);

  Face_handle intersect_arcs(Halfedge_handle ha, Halfedge_handle hb);
  Face_handle unintersect_arcs(Halfedge_handle ha, Halfedge_handle hb);


  void delete_component(Vertex_handle vh) ;

  void delete_circle(Vertex_handle v) ;

  void audit(bool extra_vertices=false) const ;	 


  typedef std::pair<Point, Curve> NFP;

  void reserve(int nv, int ne, int nf);

 
  void set_has_boundary(bool tf);

  bool has_boundary() const;

  void set_curve(Halfedge_handle h, Curve c) const {
    h->set_curve(c);
    h->opposite()->set_curve(c.opposite());
  }

  Halfedge_handle remove_vertex(Vertex_handle v);

  //Halfedge_handle remove_vertex(Vertex_handle v, Curve c);

  void move_target(Halfedge_handle h, Vertex_handle nt, bool cleanup);

  void exchange_sphere_extremums(Curve::Key k, Curve::Key l);

  Halfedge_handle halfedge(Vertex_handle v, Face_handle f) const ;




  template <class HH>
  HH next_edge_on_circle(HH h) const {
    CGAL_precondition(h== HH() || is_valid(h));
    CGAL_precondition(h->curve().is_arc());
    HH r= h->next();
    do {
      if (h->curve().is_arc() == r->curve().is_arc()
	  && h->curve().key() == r->curve().key()) {
	CGAL_assertion(h->curve().is_inside()
		     == r->curve().is_inside());
	return r;
      }
      r=r->opposite()->next();
      CGAL_assertion(r!= h->opposite());
    } while (true);
    CGAL_error();
    return HH();
  }




  template <class HH>
  HH next_edge_on_rule(HH h) const {
    CGAL_precondition(h== HH() || is_valid(h));
    CGAL_precondition(h->curve().is_rule());
    if (!h->vertex()->point().is_rule_rule()) return HH();
    //int deg = degree(h->vertex());
    HH r= h->next();
    do {
      if (r->curve().constant_coordinate() == h->curve().constant_coordinate()) {
	return r;
      }
      r=r->opposite()->next();
    } while (r != h->opposite());
    return HH();
  }

  // an outward facing edge
  Halfedge_handle cross_edge(Halfedge_handle) const;

  // a halfedge on the curve (an inside one)
  Halfedge_handle a_halfedge(Curve::Key k) const;

  // a halfedge on the rule pointing to the extremal vertex
  // Halfedge_handle rule_halfedge(Curve::Key k, Rule_direction i) const;

  // a halfedge on the circle pointing to the extremal vertex (inside)
  //Halfedge_handle extremum_halfedge(Curve::Key k, Rule_direction i) const;
  
  // insert the vertex so that h->opposite points to it
  //Halfedge_handle insert_vertex( Point p, Halfedge_handle h);


  /*std::pair<Halfedge_handle,Halfedge_handle> 
  intersect(Halfedge_handle ha, Halfedge_handle hb);

  std::pair<Halfedge_handle,Halfedge_handle>  unintersect(Face_handle f);*/

 
  
  //friend class Cross_section_initializer CGAL_AOS3_TARG;
 

  //void relabel_rule(Halfedge_handle h, Curve nl);

 

  //typedef CGAL::array<Halfedge_handle, 4> Halfedge_quad;

  /*const Halfedge_handle halfedges(CGAL_AOS3_TYPENAME Curve::Key k) const {
    return halfedges_[k.input_index()];
    }

    Halfedge_quad& halfedges(CGAL_AOS3_TYPENAME Curve::Key k)  {
    return halfedges_[k.input_index()];
    }*/
  void set_number_of_spheres(unsigned int i) {
    CGAL_assertion(halfedges_.size() <=i);
    halfedges_.resize(i);
  }

  void new_target(CGAL_AOS3_TYPENAME Curve::Key k, Halfedge_handle c[]);

  Vertex_handle new_vertex(Point p);
  

  // insert the vertex so that h->opposite points to it
  Halfedge_handle insert_vertex(Vertex_handle  vh, Halfedge_handle h);

  bool is_valid(Halfedge_const_handle h) const {
    for (Halfedge_const_iterator it =halfedges_begin(); it != halfedges_end(); ++it) {
      if (it==h) return true;
    }
    return false;
  }

  bool is_valid(Vertex_const_handle h) const {
    for (Vertex_const_iterator it =vertices_begin(); it != vertices_end(); ++it) {
      if (it==h) return true;
    }
    return false;
  }

  bool is_valid(Face_const_handle h) const {
    for (Face_const_iterator it =faces_begin(); it != faces_end(); ++it) {
      if (it==h) return true;
    }
    return false;
  }

  void swap_labels(Curve::Key k, Curve::Key l);

private:

  void swap_curves(Halfedge_handle k, Halfedge_handle l);
  

  void connect(Halfedge_handle a, Halfedge_handle b);

  /*typedef boost::tuple<Halfedge_handle, Halfedge_handle,
    Halfedge_handle, Halfedge_handle> Halfedge_handle_quadruple;*/

  //void new_target(Curve::Key k, Halfedge_handle tar[]);
  
  void set_halfedge(Halfedge_handle h);


  /*Halfedge_handle split_face(Halfedge_handle o, Halfedge_handle d,
    Curve c);*/

  //void relabel_target(Halfedge_handle v[], Curve::Key k);

  
  //void exchange_vertices(Halfedge_handle h, Halfedge_handle p);

  void audit_vertex(Vertex_const_handle v, bool extra) const ;
  


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
  
  


  // insert a vertex for p in both edges
  // merge the two vertices
  /*std::pair<Halfedge_handle, Halfedge_handle> pinch_bl(Halfedge_handle a,
						       Halfedge_handle b, Point p);

  // a and b go to the same vertex. Make them go to different vertices
  // return the new vertex on b
  Vertex_handle unpinch_bl(Halfedge_handle a, Halfedge_handle b);

  void merge_vertices_bl(Halfedge_handle a, Halfedge_handle b);*/

  Halfedge_handle new_halfedge(Curve c);

  //typedef std::pair<Point,Point>  ED;

  
 

 
  mutable HDS hds_;
  Face_handle inf_;
  mutable std::vector<Curve> errors_;
  //std::vector<Vertex_handle> targets_;
  std::vector<Halfedge_handle> halfedges_;
  int num_components_;
  Visitor v_;

  // for new_face when constructing things
  //std::map<Edge, Halfedge_handle> unmatched_hedges_;
  //std::map<Point, Vertex_handle> points_;
    
};

CGAL_AOS3_END_INTERNAL_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Combinatorial_cross_section_impl.h>
#endif

#endif
