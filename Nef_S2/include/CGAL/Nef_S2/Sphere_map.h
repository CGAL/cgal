// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_SPHERE_MAP_H
#define CGAL_SPHERE_MAP_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_2/Object_handle.h>
#include <CGAL/Nef_S2/SM_items.h>
#include <CGAL/Nef_S2/SM_list.h>
#include <CGAL/Nef_S2/SM_iteration.h>
#include <CGAL/Nef_S2/Generic_handle_map.h>
#include <CGAL/Nef_2/iterator_tools.h>
#include <list>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 109
#include <CGAL/Nef_2/debug.h>

#include<boost/optional.hpp>
#include<boost/none.hpp>

namespace CGAL {

template <typename HE>
class move_edge_around_svertex {
public:
  void forward(HE& e) const  { e = (e->sprev()->twin()); }
  void backward(HE& e) const { e = (e->twin()->snext()); }
};

template <typename HE>
struct move_edge_around_sface {
  void forward(HE& e)  const { e = (e->snext()); }
  void backward(HE& e) const { e = (e->sprev()); }
};

/*{\Manpage {Sphere_map}{Kernel}{Sphere Maps}{M}}*/

template <typename Kernel_, typename Items_, typename Mark_>
class Sphere_map {

/*{\Mdefinition selective sphere map container based on
the HDS design of Kettner.}*/

public:
  /*{\Mtypes 7}*/
  typedef Sphere_map<Kernel_, Items_,Mark_>   Self;
  typedef Kernel_                             Sphere_kernel;
  typedef Items_                              Items;
  typedef Mark_                               Mark;

  friend class SM_const_decorator<Self>;
  friend class SM_decorator<Self>;

  typedef typename Sphere_kernel::Sphere_point     Sphere_point;
  /*{\Mtypemember points on the unit sphere.}*/
  typedef typename Sphere_kernel::Sphere_segment   Sphere_segment;
  /*{\Mtypemember segments on the unit sphere.}*/
  typedef typename Sphere_kernel::Sphere_circle    Sphere_circle;
  /*{\Mtypemember segments on the unit sphere.}*/
  typedef typename Sphere_kernel::Sphere_direction Sphere_direction;
  /*{\Mtypemember directions on the unit sphere.}*/
  //  typedef bool                                     Mark;
  /*{\Mtypemember selective attributes of all objects.}*/
  typedef size_t   Size_type;
  /*{\Mtypemember size type.}*/

  /*{\Mtext For all objects |Vertex|, |Halfedge|, |Halfloop|, |Face|
  there are handle and iterator types |xxx_handle|, |xxx_iterator|.
  There's no type |SLoop_iterator|, as there is
  at most one |SLoop| pair per sphere map.}*/

  typedef typename Items::template SVertex<Self>        SVertex_base;
  typedef SNC_in_place_list_svertex<SVertex_base>       SVertex;
  typedef CGAL::In_place_list<SVertex,false>            SVertex_list;
  typedef CGAL_ALLOCATOR(SVertex)                       SVertex_alloc;
  typedef typename SVertex_list::iterator               SVertex_handle;
  typedef typename SVertex_list::const_iterator         SVertex_const_handle;
  typedef typename SVertex_list::iterator               SVertex_iterator;
  typedef typename SVertex_list::const_iterator         SVertex_const_iterator;

  typedef typename Items::template SHalfedge<Self>      SHalfedge_base;
  typedef SNC_in_place_list_shalfedge<SHalfedge_base>   SHalfedge;
  typedef CGAL::In_place_list<SHalfedge,false>          SHalfedge_list;
  typedef CGAL_ALLOCATOR(SHalfedge)                     SHalfedge_alloc;
  typedef typename SHalfedge_list::iterator             SHalfedge_handle;
  typedef typename SHalfedge_list::const_iterator       SHalfedge_const_handle;
  typedef typename SHalfedge_list::iterator             SHalfedge_iterator;
  typedef typename SHalfedge_list::const_iterator       SHalfedge_const_iterator;

  typedef typename Items::template SFace<Self>          SFace_base;
  typedef SNC_in_place_list_sface<SFace_base>           SFace;
  typedef CGAL::In_place_list<SFace,false>              SFace_list;
  typedef CGAL_ALLOCATOR(SFace)                         SFace_alloc;
  typedef typename SFace_list::iterator                 SFace_handle;
  typedef typename SFace_list::const_iterator           SFace_const_handle;
  typedef typename SFace_list::iterator                 SFace_iterator;
  typedef typename SFace_list::const_iterator           SFace_const_iterator;

  typedef typename Items::template SHalfloop<Self>      SHalfloop;
  typedef SHalfloop*                                    SHalfloop_handle;
  typedef const SHalfloop*                              SHalfloop_const_handle;
  typedef SHalfloop*                                    SHalfloop_iterator;
  typedef const SHalfloop*                              SHalfloop_const_iterator;

  typedef CGAL::Object_handle Object_handle;
  /*{\Mtypemember a generic handle to an object of |\Mvar|. 
  The kind of the object can be determined and the object assigned 
  by the function:\\ 
  |bool assign(xxx_handle& h, Object_handle o)|\\ 
  where the function returns |true| iff the assignment of |o| to 
  |h| was valid.}*/

  typedef std::list<Object_handle>                     Object_list;
  typedef typename Object_list::iterator               Object_iterator;
  typedef typename Object_list::const_iterator         Object_const_iterator;
  typedef boost::optional<Object_iterator>             Optional_object_iterator ;
  typedef Generic_handle_map<Optional_object_iterator> Handle_to_iterator_map;

  typedef Sphere_map*       Constructor_parameter;
  typedef const Sphere_map* Constructor_const_parameter;

  class SFace_cycle_iterator : public Object_iterator 
  /*{\Mtypemember a generic iterator to an object in the boundary
  of a facet. Convertible to |Object_handle|.}*/
  { typedef Object_iterator Ibase;
  public:
    SFace_cycle_iterator() : Ibase() {}
    SFace_cycle_iterator(const Ibase& b) : Ibase(b) {}
    SFace_cycle_iterator(const SFace_cycle_iterator& i) : Ibase(i) {}  
    bool is_svertex() const 
    { SVertex_handle v; return CGAL::assign(v,Ibase::operator*()); }
    bool is_shalfedge() const
    { SHalfedge_handle e; return CGAL::assign(e,Ibase::operator*()); }
    bool is_shalfloop() const
    { SHalfloop_handle l; return CGAL::assign(l,Ibase::operator*()); }
    operator SVertex_handle() const 
    { SVertex_handle v; CGAL::assign(v,Ibase::operator*()); return v; }
    operator SHalfedge_handle() const 
    { SHalfedge_handle e; CGAL::assign(e,Ibase::operator*()); return e; }
    operator SHalfloop_handle() const 
    { SHalfloop_handle l; CGAL::assign(l,Ibase::operator*()); return l; }

    operator Object_handle() const { return Ibase::operator*(); }
    Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const 
    { CGAL_error_msg("not impl."); return Object_handle(); }
  };

  class SFace_cycle_const_iterator : public Object_const_iterator 
  /*{\Mtypemember a generic iterator to an object in the boundary
  of a facet. Convertible to |Object_handle|.}*/
  { typedef Object_const_iterator Ibase;
  public:
    SFace_cycle_const_iterator() : Ibase() {}
    SFace_cycle_const_iterator(const Ibase& b) : Ibase(b) {}
    SFace_cycle_const_iterator(const SFace_cycle_const_iterator& i) 
      : Ibase(i) {}  
    bool is_svertex() const 
    { SVertex_handle v; return CGAL::assign(v,Ibase::operator*()); }
    bool is_shalfedge() const
    { SHalfedge_handle e; return CGAL::assign(e,Ibase::operator*()); }
    bool is_shalfloop() const
    { SHalfloop_handle l; return CGAL::assign(l,Ibase::operator*()); }
    operator SVertex_const_handle() const 
    { SVertex_handle v; CGAL::assign(v,Ibase::operator*()); 
      return SVertex_const_handle(v); }
    operator SHalfedge_const_handle() const 
    { SHalfedge_handle e; CGAL::assign(e,Ibase::operator*()); 
      return SHalfedge_const_handle(e); }
    operator SHalfloop_const_handle() const 
    { SHalfloop_handle l; CGAL::assign(l,Ibase::operator*()); 
      return SHalfloop_const_handle(l); }

    operator Object_handle() const { return Ibase::operator*(); }
    const Object_handle& operator*() const { return Ibase::operator*(); }
    Object_handle  operator->() const 
    { CGAL_error_msg("not impl."); return Object_handle(); }
  };

  /*{\Mtext Local types are handles, iterators and circulators of the
    following kind: |SVertex_handle|, |SVertex_iterator|, |SHalfedge_handle|,
    |SHalfedge_iterator|, |SHalfloop_handle|, |SHalfloop_iterator|,
    |SFace_handle|, |SFace_iterator|.  Additionally the following
    circulators are defined.}*/

  typedef CircFromIt<
    SHalfedge_const_iterator, 
    move_edge_around_svertex<SHalfedge_const_iterator> > 
    SHalfedge_around_svertex_const_circulator;
  /*{\Mtypemember circulating the adjacency list of an vertex |v|.}*/
  
  typedef CircFromIt<
    SHalfedge_const_iterator, 
    move_edge_around_sface<SHalfedge_const_iterator> > 
    SHalfedge_around_sface_const_circulator;
  /*{\Mtypemember circulating the face cycle of an face |f|.}*/

  typedef CircFromIt<
    SHalfedge_iterator, 
    move_edge_around_svertex<SHalfedge_iterator> > 
    SHalfedge_around_svertex_circulator;
  /*{\Mtypemember circulating the adjacency list of an vertex |v|.}*/
  
  typedef CircFromIt<
    SHalfedge_iterator, 
    move_edge_around_sface<SHalfedge_iterator> > 
    SHalfedge_around_sface_circulator;
  /*{\Mtypemember circulating the face cycle of an face |f|.}*/
  
  /*{\Mcreation 3}*/
  /*{\Mtext |\Mname| is default and copy constructible. Note that copy
  construction means cloning an isomorphic structure and is thus an
  expensive operation.}*/

  Sphere_map(bool = false) : boundary_item_(boost::none), 
    svertices_(), sedges_(), sfaces_(), shalfloop_() {}

  ~Sphere_map() { clear(); }

  Sphere_map(const Self& D) : boundary_item_(boost::none),
    svertices_(D.svertices_), 
    sedges_(D.sedges_), 
    sfaces_(D.sfaces_), 
    shalfloop_(0)
  { if ( D.shalfloop_ != 0 ) new_shalfloop_pair(*(D.shalfloop_));
    pointer_update(D); 
  }

  Self& operator=(const Self& D) 
  { if ( this == &D ) return *this; 
    clear();
    svertices_ = D.svertices_; 
    sfaces_ = D.sfaces_;
    sedges_ = D.sedges_;
    if ( D.shalfloop_ != 0 ) new_shalfloop_pair(*D.shalfloop_);
    pointer_update(D);
    return *this;
  }

  void clear()
  { 
    boundary_item_.clear(boost::none);
    svertices_.destroy(); 
    sfaces_.destroy();
    while ( shalfedges_begin() != shalfedges_end() )
      delete_shalfedge_pair( shalfedges_begin() );
    if ( shalfloop_ != 0 ) { delete_shalfloop_pair(); shalfloop_=0; }
  }

  template <typename H>
  bool is_sm_boundary_object(H h) const
  { return boundary_item_[h]!=boost::none; }

  template <typename H>
  Object_iterator& sm_boundary_item(H h)
  { return *boundary_item_[h]; }

  template <typename H>
  void store_sm_boundary_item(H h, Object_iterator o)
  { boundary_item_[h] = o; }

  template <typename H>
  void undef_sm_boundary_item(H h)
  { CGAL_assertion(boundary_item_[h]!=boost::none);
    boundary_item_[h] = boost::none; }

  void reset_sm_iterator_hash(Object_iterator it)
  { SVertex_handle sv;
    SHalfedge_handle se;
    SHalfloop_handle sl;
    if ( CGAL::assign(se,*it) ) { undef_sm_boundary_item(se); return; }
    if ( CGAL::assign(sl,*it) ) { undef_sm_boundary_item(sl); return; }
    if ( CGAL::assign(sv,*it) ) { undef_sm_boundary_item(sv); return; }
  }

  void reset_sm_object_list(Object_list& L)
  { Object_iterator oit;
    CGAL_forall_iterators(oit,L) reset_sm_iterator_hash(oit);
    L.clear();
  }

  /*{\Moperations 2.5 3}*/

  // The constant iterators and circulators.
  SVertex_const_iterator   svertices_begin()  const { return svertices_.begin();}
  SVertex_const_iterator   svertices_end()    const { return svertices_.end();}
  SHalfedge_const_iterator shalfedges_begin() const { return sedges_.begin();}
  SHalfedge_const_iterator shalfedges_end()   const { return sedges_.end();}
  SHalfloop_const_iterator shalfloops_begin() const { return shalfloop_; }
  SHalfloop_const_iterator shalfloops_end()   const 
    { return shalfloop_ != 0 ? shalfloop_+2 : shalfloop_; }
  SFace_const_iterator     sfaces_begin()     const { return sfaces_.begin();}
  SFace_const_iterator     sfaces_end()       const { return sfaces_.end();}

  SVertex_iterator   svertices_begin()   { return svertices_.begin();}
  SVertex_iterator   svertices_end()     { return svertices_.end();}
  SHalfedge_iterator shalfedges_begin()  { return sedges_.begin();}
  SHalfedge_iterator shalfedges_end()    { return sedges_.end();}
  SHalfloop_iterator shalfloops_begin()  { return shalfloop_; }
  SHalfloop_iterator shalfloops_end()  
    { return shalfloop_ != 0 ? shalfloop_+2 : shalfloop_; }
  SFace_iterator     sfaces_begin()      { return sfaces_.begin();}
  SFace_iterator     sfaces_end()        { return sfaces_.end();}

  SFace_cycle_const_iterator sface_cycles_begin(SFace_const_handle f) const
  { return f->sface_cycles_begin(); }
  SFace_cycle_const_iterator sface_cycles_end(SFace_const_handle f) const
  { return f->sface_cycles_end(); }
  SFace_cycle_iterator sface_cycles_begin(SFace_handle f) const
  { return f->sface_cycles_begin(); }
  SFace_cycle_iterator sface_cycles_end(SFace_handle f) const
  { return f->sface_cycles_end(); }

  /*{\Mtext The list of all objects can be accessed via iterator ranges.
  For comfortable iteration we also provide iterations macros. 
  The iterator range access operations are of the following kind:\\
  |SVertex_iterator vertices_begin()/vertices_end()|\\
  |SHalfedge_iterator halfedges_begin()/halfedges_end()|\\
  |SHalfloop_iterator halfloops_begin()/halfloops_end()|\\
  |SFace_iterator faces_begin()/faces_end()| */

  Size_type number_of_svertices() const  { return svertices_.size();}
  /*{\Mop returns the number of vertices.}*/
  Size_type number_of_shalfedges() const { return sedges_.size();}
  /*{\Mop returns the number of (directed edges).}*/
  Size_type number_of_sfaces() const     { return sfaces_.size();}
  /*{\Mop returns the number of facets.}*/
  Size_type number_of_shalfloops() const 
  { return shalfloop_!=SHalfloop_handle() ? 2 : 0; }
  /*{\Mop returns the number of shalfloop.}*/

  bool empty() const 
  { return number_of_svertices() == 0 &&
      number_of_shalfedges() == 0 &&
      number_of_shalfloops() == 0 &&
      number_of_sfaces() == 0;
  }

  bool has_shalfloop() const {
    return shalfloop_ != 0;
  }

  SHalfloop_handle& shalfloop() { 
    return shalfloop_; 
  }

  SHalfloop_const_handle shalfloop() const { 
    return shalfloop_; 
  }

  SVertex_alloc vertex_allocator;
  SVertex* get_vertex_node( const SVertex& ) {
    SVertex* p = vertex_allocator.allocate(1);
    vertex_allocator.construct( p, SVertex());
    return p;
  }
  void put_vertex_node( SVertex* p) {
    vertex_allocator.destroy(p);
    vertex_allocator.deallocate( p, 1);
  }

  SHalfedge_alloc halfedge_allocator;
  SHalfedge* get_halfedge_node( const SHalfedge& ) {
    SHalfedge* p = halfedge_allocator.allocate(1);
    halfedge_allocator.construct( p, SHalfedge());
    return p;
  }
  void put_halfedge_node( SHalfedge* p) {
    halfedge_allocator.destroy(p);
    halfedge_allocator.deallocate( p, 1);
  }

  SFace_alloc face_allocator;
  SFace* get_face_node( const SFace& ) {
    SFace* p = face_allocator.allocate(1);
    face_allocator.construct( p, SFace());
    return p;
  }
  void put_face_node( SFace* p) {
    face_allocator.destroy(p);
    face_allocator.deallocate( p, 1);
  }

  SVertex_handle new_svertex(const Sphere_point& p, 
			   Mark m = Mark())
  /*{\Mop returns a new vertex at point |p| marked by |m|.}*/
  { SVertex_handle vh = new_svertex(); vh->point() = p; vh->mark() = m;
    CGAL_NEF_TRACEN("new_svertex "<<&*vh);
    return vh;
  }

  template <typename H>
  void make_twins(H h1, H h2) { h1->twin() = h2; h2->twin() = h1; }

  SVertex_handle new_svertex() { 
    svertices_.push_back( * get_vertex_node(SVertex())); 
    return --svertices_end(); 
  }

  SFace_handle new_sface() { 
    sfaces_.push_back( * get_face_node(SFace())); 
    return --sfaces_end(); 
  } 

  SHalfedge_handle new_shalfedge_pair() { 
    SHalfedge* ep2 = get_halfedge_node(SHalfedge());
    SHalfedge* ep1 = get_halfedge_node(SHalfedge());
    sedges_.push_back( *ep1 );
    SHalfedge_handle e1 = --shalfedges_end();
    sedges_.push_back( *ep2 );
    SHalfedge_handle e2 = --shalfedges_end();
    make_twins(e1,e2); return e1; }

  SHalfloop_handle new_shalfloop_pair()
  { SHalfloop_handle ph = new SHalfloop[2];
    SHalfloop* pt(ph); ++pt;
    make_twins(ph,pt);
    shalfloop_=ph; return ph; }

  SHalfedge_handle new_shalfedge_pair(const SHalfedge& e1)
  { const SHalfedge& e2 = *(e1.twin());
    SHalfedge* ep2 = new SHalfedge[2];
    SHalfedge* ep1 = ep2++;
    *ep1=e1; *ep2=e2;
    sedges_.push_back( *ep1 );
    SHalfedge_handle eh1 = --shalfedges_end();
    sedges_.push_back( *ep2 );
    SHalfedge_handle eh2 = --shalfedges_end();
    make_twins(eh1,eh2); return eh1; }

  SHalfloop_handle new_shalfloop_pair(const SHalfloop& l1)
  { const SHalfloop& l2 = *(l1.twin());
    SHalfloop* ph = new SHalfloop[2];
    SHalfloop* pt(ph); ++pt;
    *ph=l1; *pt=l2; make_twins(ph,pt);
    shalfloop_=ph; return ph; }

  void delete_svertex(SVertex_handle h) { 
    svertices_.erase(h); 
    put_vertex_node(&*h); 
  }

  void delete_sface(SFace_handle h) { 
    sfaces_.erase(h); 
    put_face_node(&*h);
  }

  void delete_shalfedge_pair(SHalfedge_handle h) { 
    SHalfedge_handle t = h->twin();
    sedges_.erase(h); sedges_.erase(t);
    put_halfedge_node(&*h);
    put_halfedge_node(&*t);
  }

  void delete_shalfloop_pair() { 
    CGAL_assertion(has_shalfloop());
    SHalfloop* ph = &*shalfloop_;
    SHalfloop* pt = &*shalfloop_->twin();
    if ( ph > pt ) std::swap(ph,pt);
    shalfloop_ = SHalfloop_handle();
    delete [] ph; }

protected:
  void pointer_update(const Self& D);
  Handle_to_iterator_map boundary_item_;

  SVertex_list       svertices_;
  SHalfedge_list     sedges_;
  SFace_list         sfaces_;
  SHalfloop_iterator shalfloop_;

}; // Sphere_map


template <typename K, typename I, typename M>
void Sphere_map<K, I, M>::
pointer_update(const Sphere_map<K, I, M>& D)
{
  CGAL::Unique_hash_map<SVertex_const_handle,SVertex_handle>     VM;
  CGAL::Unique_hash_map<SHalfedge_const_handle,SHalfedge_handle> EM;
  CGAL::Unique_hash_map<SHalfloop_const_handle,SHalfloop_handle> LM;
  CGAL::Unique_hash_map<SFace_const_handle,SFace_handle>         FM;

  SVertex_const_iterator vc = D.svertices_begin();
  SVertex_iterator v = svertices_begin();
  for ( ; vc != D.svertices_end(); ++vc,++v) VM[vc] = v;
  VM[D.svertices_end()] = svertices_end();

  SHalfedge_const_iterator ec = D.shalfedges_begin();
  SHalfedge_iterator e = shalfedges_begin();
  for ( ; ec != D.shalfedges_end(); ++ec,++e) { 
    EM[ec] = e;
    e->mark() = ec->mark();
  }
  EM[D.shalfedges_end()] = shalfedges_end();

  SFace_const_iterator fc = D.sfaces_begin();
  SFace_iterator f = sfaces_begin();
  for ( ; fc != D.sfaces_end(); ++fc,++f) FM[fc] = f;
  FM[D.sfaces_end()] = sfaces_end();

  SHalfloop_iterator l;
  if ( D.shalfloop_ != 0 ) {
    LM[D.shalfloop_] = shalfloop_;
    LM[D.shalfloop_->twin()] = shalfloop_->twin();
    l = shalfloop();
    shalfloop_->mark() = D.shalfloop()->mark();
    shalfloop_->twin()->mark() = D.shalfloop()->twin()->mark();    
    if( !l->is_twin() && D.shalfloop()->is_twin()) l->mark() = l->twin()->mark();
  }

  for (v = svertices_begin(); v != svertices_end(); ++v) {
    // Local Graph update: (SVertices are postponed/updated as Edges)
    v->out_sedge() = EM[v->out_sedge()];
    v->incident_sface() = FM[v->incident_sface()];
  }
  // Edge update:
  for (e = shalfedges_begin(); e != shalfedges_end(); ++e) {
    e->twin() = EM[e->twin()];
    e->sprev() = EM[e->sprev()];
    e->snext() = EM[e->snext()];
    e->source() = VM[e->source()];
    e->incident_sface() = FM[e->incident_sface()];
  }
  for ( ec = D.shalfedges_begin(), e = shalfedges_begin();
	ec != D.shalfedges_end(); ++ec, ++e) {
    if( !e->is_twin() && ec->is_twin()) e->mark() = e->twin()->mark();
  }

  for (l = shalfloops_begin(); l != shalfloops_end(); ++l) {
    //    l->twin() = LM[l->twin()];
    l->incident_sface() = FM[l->incident_sface()];
  }


  for (f = sfaces_begin(); f != sfaces_end(); ++f) {
    SFace_cycle_iterator fci; 
    for(fci = f->boundary_entry_objects().begin();
	fci != f->boundary_entry_objects().end(); ++fci) {
      if ( fci.is_svertex() ) 
      { v = SVertex_handle(fci);
	*fci = make_object(VM[v]); store_sm_boundary_item(v,fci); }
      else if ( fci.is_shalfedge() ) 
      { e = SHalfedge_handle(fci);
	*fci = make_object(EM[e]); store_sm_boundary_item(e,fci); }
      else if ( fci.is_shalfloop() ) 
      { l = SHalfloop_handle(fci);
	*fci = make_object(LM[l]); store_sm_boundary_item(l,fci); }
      else CGAL_error_msg("damn wrong boundary item in face.");
    }
  }
}


} //namespace CGAL
#endif // CGAL_SPHERE_MAP_H
