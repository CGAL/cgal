// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_S2/Sphere_map.h
// package       : Nef_S2 
// chapter       : Nef Polyhedra on the sphere
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
// ============================================================================

#ifndef CGAL_SPHERE_MAP_H
#define CGAL_SPHERE_MAP_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_S2/nef_assertions.h>
#include <CGAL/Nef_S2/SM_items.h>
#include <CGAL/Nef_S2/SM_iteration.h>
#include <CGAL/Nef_S2/Generic_handle_map.h>
#include <list>
#undef _DEBUG
#define _DEBUG 109
#include <CGAL/Nef_S2/debug.h>

CGAL_BEGIN_NAMESPACE

/*{\Manpage {Sphere_map}{Kernel}{Sphere Maps}{M}}*/

template <typename Kernel_>
class Sphere_map {

/*{\Mdefinition selective sphere map container based on
the HDS design of Kettner.}*/

public:
  /*{\Mtypes 7}*/
  typedef Kernel_ Kernel;
  typedef Sphere_map<Kernel_>    Self;
  typedef SM_items<Kernel_,bool> Items;

  friend class SM_const_decorator<Self,Kernel>;
  friend class SM_decorator<Self,Kernel>;

  typedef typename Kernel::Sphere_point Sphere_point;
  /*{\Mtypemember embedding vertices.}*/
  typedef typename Kernel::Sphere_circle Sphere_circle;
  /*{\Mtypemember embedding edges.}*/
  typedef bool Mark;
  /*{\Mtypemember selective attributes of all objects.}*/
  typedef size_t   Size_type;
  /*{\Mtypemember size type.}*/

  /*{\Mtext For all objects |Vertex|, |Halfedge|, |Halfloop|, |Face|
  there are handle and iterator types |xxx_handle|, |xxx_iterator|.
  There's no type |SLoop_iterator|, as there is
  at most one |SLoop| pair per sphere map.}*/

  typedef typename Items::template Vertex<Self>    Vertex;
  typedef CGAL::In_place_list<Vertex,false>        Vertex_list;
  typedef typename Vertex_list::iterator           Vertex_handle;
  typedef typename Vertex_list::const_iterator     Vertex_const_handle;
  typedef typename Vertex_list::iterator           Vertex_iterator;
  typedef typename Vertex_list::const_iterator     Vertex_const_iterator;

  typedef typename Items::template Halfedge<Self>  Halfedge;
  typedef CGAL::In_place_list<Halfedge,false>      Halfedge_list;
  typedef typename Halfedge_list::iterator         Halfedge_handle;
  typedef typename Halfedge_list::const_iterator   Halfedge_const_handle;
  typedef typename Halfedge_list::iterator         Halfedge_iterator;
  typedef typename Halfedge_list::const_iterator   Halfedge_const_iterator;

  typedef typename Items::template Face<Self>      Face;
  typedef CGAL::In_place_list<Face,false>          Face_list;
  typedef typename Face_list::iterator             Face_handle;
  typedef typename Face_list::const_iterator       Face_const_handle;
  typedef typename Face_list::iterator             Face_iterator;
  typedef typename Face_list::const_iterator       Face_const_iterator;

  typedef typename Items::template Halfloop<Self>  Halfloop;
  typedef Halfloop* Halfloop_handle;
  typedef const Halfloop* Halfloop_const_handle;
  typedef Halfloop* Halfloop_iterator;
  typedef const Halfloop* Halfloop_const_iterator;

  class Object_handle 
  /*{\Mtypemember a generic handle to an object of |\Mvar|. 
  The kind of the object can be determined and the object assigned 
  by the function:\\ 
  |bool assign(xxx_handle& h, Object_handle o)|\\ 
  where the function returns |true| iff the assignment of |o| to 
  |h| was valid.}*/
    : public CGAL::Object 
  {
    typedef CGAL::Object Base;
  public:
    Object_handle() : Base() {}
    Object_handle(const CGAL::Object& o) : Base(o) {}
    Object_handle(const Object_handle& h) : Base(h) {}
    Object_handle& operator=(const Object_handle& h) 
    { Base::operator=(h); return *this; }
    bool operator==(CGAL_NULL_TYPE n) const
    { assert(n == 0); return Base::is_empty(); }
    bool operator!=(CGAL_NULL_TYPE n) const
    { assert(n == 0); return !Base::is_empty(); }
  }; // Object_handle

  typedef std::list<Object_handle>    Object_list;
  typedef typename Object_list::iterator       Object_iterator;
  typedef typename Object_list::const_iterator Object_const_iterator;
  typedef Generic_handle_map<Object_iterator> Handle_to_iterator_map;

  class Face_cycle_iterator : public Object_iterator 
  /*{\Mtypemember a generic iterator to an object in the boundary
  of a facet. Convertible to |Object_handle|.}*/
  { typedef Object_iterator Ibase;
  public:
    Face_cycle_iterator() : Ibase() {}
    Face_cycle_iterator(const Ibase& b) : Ibase(b) {}
    Face_cycle_iterator(const Face_cycle_iterator& i) : Ibase(i) {}  
    bool is_vertex() const 
    { Vertex_handle v; return CGAL::assign(v,Ibase::operator*()); }
    bool is_halfedge() const
    { Halfedge_handle e; return CGAL::assign(e,Ibase::operator*()); }
    bool is_halfloop() const
    { Halfloop_handle l; return CGAL::assign(l,Ibase::operator*()); }
    operator Vertex_handle() const 
    { Vertex_handle v; CGAL::assign(v,Ibase::operator*()); return v; }
    operator Halfedge_handle() const 
    { Halfedge_handle e; CGAL::assign(e,Ibase::operator*()); return e; }
    operator Halfloop_handle() const 
    { Halfloop_handle l; CGAL::assign(l,Ibase::operator*()); return l; }

    operator Object_handle() const { return Ibase::operator*(); }
    Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const { CGAL_nef_assertion_msg(0,"not impl."); }
  };

  class Face_cycle_const_iterator : public Object_const_iterator 
  /*{\Mtypemember a generic iterator to an object in the boundary
  of a facet. Convertible to |Object_handle|.}*/
  { typedef Object_const_iterator Ibase;
  public:
    Face_cycle_const_iterator() : Ibase() {}
    Face_cycle_const_iterator(const Ibase& b) : Ibase(b) {}
    Face_cycle_const_iterator(const Face_cycle_const_iterator& i) 
      : Ibase(i) {}  
    bool is_vertex() const 
    { Vertex_handle v; return CGAL::assign(v,Ibase::operator*()); }
    bool is_halfedge() const
    { Halfedge_handle e; return CGAL::assign(e,Ibase::operator*()); }
    bool is_halfloop() const
    { Halfloop_handle l; return CGAL::assign(l,Ibase::operator*()); }
    operator Vertex_const_handle() const 
    { Vertex_handle v; CGAL::assign(v,Ibase::operator*()); 
      return Vertex_const_handle(v); }
    operator Halfedge_const_handle() const 
    { Halfedge_handle e; CGAL::assign(e,Ibase::operator*()); 
      return Halfedge_const_handle(e); }
    operator Halfloop_const_handle() const 
    { Halfloop_handle l; CGAL::assign(l,Ibase::operator*()); 
      return Halfloop_const_handle(l); }

    operator Object_handle() const { return Ibase::operator*(); }
    const Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const { CGAL_nef_assertion_msg(0,"not impl."); }
  };

  /*{\Mcreation 3}*/
  /*{\Mtext |\Mname| is default and copy constructible. Note that copy
  construction means cloning an isomorphic structure and is thus an
  expensive operation.}*/

  Sphere_map() : boundary_item_(undef_), 
    vertices_(), edges_(), faces_(), loops_(0) 
  { m_pos_ = m_neg_ = Mark(); }

  ~Sphere_map() { clear(); }

  Sphere_map(const Self& D) : boundary_item_(undef_),
    vertices_(D.vertices_), 
    // edges_(D.edges_), 
    faces_(D.faces_), 
    loops_(0)
  { if ( D.loops_ != 0 ) new_halfloop_pair(*(D.loops_));
    Halfedge_const_iterator e;
    CGAL_forall_edges(e,D) new_halfedge_pair(*e);
    pointer_update(D); 
    m_pos_ = D.m_pos_; m_neg_ = D.m_neg_; }

  Self& operator=(const Self& D) 
  { if ( this == &D ) return *this; 
    clear();
    vertices_ = D.vertices_; 
    faces_ = D.faces_;
    Halfedge_const_iterator e;
    CGAL_forall_edges(e,D) new_halfedge_pair(*e);
    if ( D.loops_ != 0 ) new_halfloop_pair(*D.loops_);
    pointer_update(D);
    m_pos_ = D.m_pos_; m_neg_ = D.m_neg_;
    return *this;
  }

  void clear()
  { 
    boundary_item_.clear(undef_);
    vertices_.destroy(); 
    faces_.destroy();
    while ( halfedges_begin() != halfedges_end() )
      delete_halfedge_pair( halfedges_begin() );
    if ( loops_ != 0 ) { delete_halfloop_pair(loops_); loops_=0; }
    m_pos_ = m_neg_ = Mark(); 
  }

  template <typename H>
  bool is_boundary_object(H h) 
  { return boundary_item_[h]!=undef_; }

  template <typename H>
  Object_iterator& boundary_item(H h)
  { return boundary_item_[h]; }

  template <typename H>
  void store_boundary_item(H h, Object_iterator o)
  { boundary_item_[h] = o; }

  template <typename H>
  void undef_boundary_item(H h)
  { CGAL_nef_assertion(boundary_item_[h]!=undef_);
    boundary_item_[h] = undef_; }

  void reset_iterator_hash(Object_iterator it)
  { Vertex_handle sv;
    Halfedge_handle se;
    Halfloop_handle sl;
    if ( assign(se,*it) ) { undef_boundary_item(se); return; }
    if ( assign(sl,*it) ) { undef_boundary_item(sl); return; }
    if ( assign(sv,*it) ) { undef_boundary_item(sv); return; }
  }

  void reset_object_list(Object_list& L)
  { Object_iterator oit;
    CGAL_forall_iterators(oit,L) reset_iterator_hash(oit);
    L.clear();
  }

  /*{\Moperations 2.5 3}*/

  // The constant iterators and circulators.
  Vertex_const_iterator   vertices_begin()  const { return vertices_.begin();}
  Vertex_const_iterator   vertices_end()    const { return vertices_.end();}
  Halfedge_const_iterator halfedges_begin() const { return edges_.begin();}
  Halfedge_const_iterator halfedges_end()   const { return edges_.end();}
  Halfloop_const_iterator halfloops_begin() const { return loops_; }
  Halfloop_const_iterator halfloops_end()   const 
    { return loops_ != 0 ? loops_+2 : loops_; }
  Face_const_iterator     faces_begin()     const { return faces_.begin();}
  Face_const_iterator     faces_end()       const { return faces_.end();}

  Vertex_iterator   vertices_begin()   { return vertices_.begin();}
  Vertex_iterator   vertices_end()     { return vertices_.end();}
  Halfedge_iterator halfedges_begin()  { return edges_.begin();}
  Halfedge_iterator halfedges_end()    { return edges_.end();}
  Halfloop_iterator halfloops_begin()  { return loops_; }
  Halfloop_iterator halfloops_end()  
    { return loops_ != 0 ? loops_+2 : loops_; }
  Face_iterator     faces_begin()      { return faces_.begin();}
  Face_iterator     faces_end()        { return faces_.end();}

  Face_cycle_const_iterator face_cycles_begin(Face_const_handle f) const
  { return f->face_cycles_begin(); }
  Face_cycle_const_iterator face_cycles_end(Face_const_handle f) const
  { return f->face_cycles_end(); }
  Face_cycle_iterator face_cycles_begin(Face_handle f) const
  { return f->face_cycles_begin(); }
  Face_cycle_iterator face_cycles_end(Face_handle f) const
  { return f->face_cycles_end(); }

  /*{\Mtext The list of all objects can be accessed via iterator ranges.
  For comfortable iteration we also provide iterations macros. 
  The iterator range access operations are of the following kind:\\
  |Vertex_iterator vertices_begin()/vertices_end()|\\
  |Halfedge_iterator halfedges_begin()/halfedges_end()|\\
  |Halfloop_iterator halfloops_begin()/halfloops_end()|\\
  |Face_iterator faces_begin()/faces_end()|

  The macros are then |CGAL_forall_vertices(v,\Mvar)|, 
  |CGAL_forall_halfedges(e,\Mvar)|,
  |CGAL_forall_halfloops(l,\Mvar)|, 
  |CGAL_forall_faces(f,\Mvar)|.}*/

  Size_type number_of_vertices() const  { return vertices_.size();}
  /*{\Mop returns the number of vertices.}*/
  Size_type number_of_halfedges() const { return edges_.size();}
  /*{\Mop returns the number of (directed edges).}*/
  Size_type number_of_faces() const     { return faces_.size();}
  /*{\Mop returns the number of facets.}*/
  Size_type number_of_halfloops() const 
  { return loops_!=Halfloop_handle() ? 2 : 0; }
  /*{\Mop returns the number of sloops.}*/

  bool empty() const 
  { return number_of_vertices() == 0 &&
      number_of_halfedges() == 0 &&
      number_of_halfloops() == 0 &&
      number_of_faces() == 0;
  }

  Vertex_handle new_vertex(const Sphere_point& p, 
			   Mark m = Mark())
  /*{\Mop returns a new vertex at point |p| marked by |m|.}*/
  { Vertex_handle vh = new_vertex(); vh->point_ = p; vh->mark_ = m;
    TRACEN("new_vertex "<<&*vh);
    return vh;
  }

  template <typename H>
  void make_twins(H h1, H h2) { h1->twin_ = h2; h2->twin_ = h1; }

  Vertex_handle new_vertex()
  { vertices_.push_back( * new Vertex); return --vertices_end(); }

  Face_handle new_face()
  { faces_.push_back( * new Face ); return --faces_end(); } 

  Halfedge_handle new_halfedge_pair()
  { Halfedge* ep2 = new Halfedge[2];
    Halfedge* ep1 = ep2++;
    edges_.push_back( *ep1 );
    Halfedge_handle e1 = --halfedges_end();
    edges_.push_back( *ep2 );
    Halfedge_handle e2 = --halfedges_end();
    make_twins(e1,e2); return e1; }

  Halfloop_handle new_halfloop_pair()
  { Halfloop_handle ph = new Halfloop[2];
    Halfloop* pt(ph); ++pt;
    make_twins(ph,pt);
    loops_=ph; return ph; }

  Halfedge_handle new_halfedge_pair(const Halfedge& e1)
  { const Halfedge& e2 = *(e1.twin_);
    Halfedge* ep2 = new Halfedge[2];
    Halfedge* ep1 = ep2++;
    *ep1=e1; *ep2=e2;
    edges_.push_back( *ep1 );
    Halfedge_handle eh1 = --halfedges_end();
    edges_.push_back( *ep2 );
    Halfedge_handle eh2 = --halfedges_end();
    make_twins(eh1,eh2); return eh1; }

  Halfloop_handle new_halfloop_pair(const Halfloop& l1)
  { const Halfloop& l2 = *(l1.twin_);
    Halfloop* ph = new Halfloop[2];
    Halfloop* pt(ph); ++pt;
    *ph=l1; *pt=l2; make_twins(ph,pt);
    loops_=ph; return ph; }

  void delete_vertex(Vertex_handle h)
  { vertices_.erase(h); delete &* h; }

  void delete_face(Face_handle h)
  { faces_.erase(h); delete &* h; }

  void delete_halfedge_pair(Halfedge_handle h)
  { Halfedge_handle t = h->twin_;
    edges_.erase(h); edges_.erase(t);
    Halfedge* ph = &*h;
    Halfedge* pt = &*t;
    if ( ph > pt ) std::swap(ph,pt);
    delete [] ph; }

  void delete_halfloop_pair(Halfloop_handle h)
  { Halfloop* ph = &*h;
    Halfloop* pt = &*(h->twin_);
    if ( ph > pt ) std::swap(ph,pt);
    loops_ = Halfloop_handle();
    delete [] ph; }

protected:
  void pointer_update(const Self& D);
  Handle_to_iterator_map boundary_item_;
  static Object_iterator undef_;

  Vertex_list       vertices_;
  Halfedge_list     edges_;
  Face_list         faces_;
  Halfloop_iterator loops_;
  Mark m_pos_, m_neg_; 
  /* two default marks at y- 
     m_neg_ just below the pos-xy-equator
     m_pos_ just above the neg-xy-equator */

}; // Sphere_map


template <typename K>
void Sphere_map<K>::
pointer_update(const Sphere_map<K>& D)
{
  CGAL::Unique_hash_map<Vertex_const_handle,Vertex_handle>     VM;
  CGAL::Unique_hash_map<Halfedge_const_handle,Halfedge_handle> EM;
  CGAL::Unique_hash_map<Halfloop_const_handle,Halfloop_handle> LM;
  CGAL::Unique_hash_map<Face_const_handle,Face_handle>         FM;

  Vertex_const_iterator vc = D.vertices_begin();
  Vertex_iterator v = vertices_begin();
  for ( ; vc != D.vertices_end(); ++vc,++v) VM[vc] = v;
  VM[D.vertices_end()] = vertices_end();

  Halfedge_const_iterator ec = D.halfedges_begin();
  Halfedge_iterator e = halfedges_begin();
  for ( ; ec != D.halfedges_end(); ++ec,++e) EM[ec] = e;
  EM[D.halfedges_end()] = halfedges_end();

  Face_const_iterator fc = D.faces_begin();
  Face_iterator f = faces_begin();
  for ( ; fc != D.faces_end(); ++fc,++f) FM[fc] = f;
  FM[D.faces_end()] = faces_end();

  Halfloop_iterator l;
  if ( D.loops_ != 0 ) {
    LM[D.loops_] = loops_;
    LM[D.loops_->twin_] = loops_->twin_;
  }

  for (v = vertices_begin(); v != vertices_end(); ++v) {
    // Local Graph update: (SVertices are postponed/updated as Edges)
    v->edge_ = EM[v->edge_];
    v->face_ = FM[v->face_];
  }
  // Edge update:
  for (e = halfedges_begin(); e != halfedges_end(); ++e) {
    // e->twin_ = EM[e->twin_]; twin is set on construction
    e->prev_ = EM[e->prev_];
    e->next_ = EM[e->next_];
    e->source_ = VM[e->source_];
    e->face_ = FM[e->face_];
  }

  for (l = halfloops_begin(); l != halfloops_end(); ++l) {
    // l->twin_ = LM[l->twin_]; twin is set on construction
    l->face_ = FM[l->face_];
  }

  for (f = faces_begin(); f != faces_end(); ++f) {
    Face_cycle_iterator fci; 
    for(fci = f->boundary_.begin(); fci != f->boundary_.end(); ++fci) {
      if ( assign(v,Object_handle(fci)) ) 
      { *fci = make_object(VM[v]); store_boundary_item(v,fci); }
      else if ( assign(e,Object_handle(fci)) ) 
      { *fci = make_object(EM[e]); store_boundary_item(e,fci); }
      else if ( assign(l,Object_handle(fci)) ) 
      { *fci = make_object(LM[l]); store_boundary_item(l,fci); }
      else CGAL_nef_assertion_msg(0,"damn wrong boundary item in face.");
    }
  }
}

template <typename Kernel_> 
typename Sphere_map<Kernel_>::Object_iterator
Sphere_map<Kernel_>::undef_;


CGAL_END_NAMESPACE
#endif // CGAL_SPHERE_MAP_H


