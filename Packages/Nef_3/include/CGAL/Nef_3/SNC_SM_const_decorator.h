// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : include/CGAL/Nef_3/SNC_SM_const_decorator.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// SNC_SM_const_decorator.h        exploration of sphere maps
// ============================================================================
#ifndef CGAL_SNC_SM_CONST_DECORATOR_H
#define CGAL_SNC_SM_CONST_DECORATOR_H

#include <CGAL/basic.h>
#include <CGAL/circulator.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_2/Object_index.h>
#include <CGAL/Nef_3/nef3_assertions.h>
#include <CGAL/Nef_3/SNC_iteration.h>
#include <string>
#include <list>
#include <sstream>
#undef _DEBUG
#define _DEBUG 67
#include <CGAL/Nef_3/debug.h>

CGAL_BEGIN_NAMESPACE

/*{\Manpage {SNC_SM_const_decorator}{Refs}
  {Topological sphere map decorator}{D}}*/

template <typename Refs_>
class SNC_SM_const_decorator 
{ typedef SNC_SM_const_decorator<Refs_> Self;
public:

/*{\Mdefinition ...}*/

/*{\Mtypes 5}*/
typedef typename Refs_::Sphere_kernel Sphere_kernel;
/*{\Mtypemember spherical geometry.}*/

typedef typename Refs_::Sphere_point   Sphere_point;
/*{\Mtypemember embedding vertices.}*/

typedef typename Refs_::Sphere_segment Sphere_segment;
/*{\Mtypemember embedding edges.}*/

typedef typename Refs_::Sphere_circle  Sphere_circle;
/*{\Mtypemember embedding loops.}*/

typedef typename Refs_::Sphere_direction Sphere_direction;
/*{\Mtypemember embedding directions.}*/

typedef typename Refs_::Mark   Mark;
/*{\Mtypemember attributes of objects (vertices, edges, faces).}*/

typedef size_t Size_type;
/*{\Mtypemember size type.}*/

typedef void*  GenPtr;

#define USING(t) typedef typename Refs_::t t
USING(Vertex_handle);
USING(Vertex_const_handle);
USING(SVertex_const_handle);
USING(SVertex_const_iterator);
USING(Vertex); 
USING(SHalfedge); 
USING(SHalfedge_const_handle); 
USING(SHalfedge_const_iterator);
USING(SHalfloop);
USING(SHalfloop_const_handle);
USING(SHalfloop_const_iterator);
USING(SFace);
USING(SFace_const_handle);
USING(SFace_const_iterator);
#undef USING

/*{\Mtext Local types are handles, iterators and circulators of the
following kind: |SVertex_handle|, |SVertex_iterator|, |SHalfedge_handle|,
|SHalfedge_iterator|, |SHalfloop_handle|, |SHalfloop_iterator|,
|SFace_handle|, |SFace_iterator|.  Additionally the following
circulators are defined.}*/

typedef typename Refs_::SHalfedge_around_svertex_const_circulator
        SHalfedge_around_svertex_const_circulator;
/*{\Mtypemember circulating the adjacency list of an vertex |v|.}*/

typedef typename Refs_::SHalfedge_around_sface_const_circulator
        SHalfedge_around_sface_const_circulator;
/*{\Mtypemember circulating the face cycle of an face |f|.}*/

typedef typename Refs_::SFace_cycle_const_iterator 
        SFace_cycle_const_iterator;
/*{\Mtypemember iterating all sface cycles of an sface |f|.
The iterator has method |bool is_svertex()|, |bool is_shalfedge()|,
|bool is_shalfloop()|, and can be converted to the corresponding
handles |SVertex_const_handle|, |SHalfedge_const_handle|, or 
|SHalfloop_const_handle|.}*/

protected: 
  Vertex_handle psm_;
  friend class SNC_SM_decorator<Refs_>;
public:

/*{\Mcreation 3}*/
SNC_SM_const_decorator() : psm_(0) {}
SNC_SM_const_decorator(const Self& D) : psm_(D.psm_) {}
Self& operator=(const Self& D) { psm_=D.psm_; return *this; }

SNC_SM_const_decorator(Vertex_handle v) : psm_(v) {}
/*{\Mcreate constructs a plane map decorator exploring the local graph
of |v|.}*/

/*{\Moperations 4 4}*/

Vertex_const_handle center_vertex() const { return psm_; }
Refs_* sncp() const { return psm_->sncp(); }

SVertex_const_handle source(SHalfedge_const_handle e) const
/*{\Mop returns the source of |e|.}*/
{ return e->source_; }

SVertex_const_handle target(SHalfedge_const_handle e) const
/*{\Mop returns the target of |e|.}*/
{ return e->twin_->source_; }

SHalfedge_const_handle twin(SHalfedge_const_handle e) const
/*{\Mop returns the twin of |e|.}*/
{ return e->twin_; }

SHalfloop_const_handle twin(SHalfloop_const_handle l) const 
/*{\Mop returns the twin of |l|.}*/
{ return l->twin_; }

bool is_isolated(SVertex_const_handle v) const
/*{\Mop returns |true| when |v| is linked to the interior of a face.}*/
{ return (SHalfedge_const_handle(v->out_sedge_) == SHalfedge_const_handle()); }

SHalfedge_const_handle first_out_edge(SVertex_const_handle v) const
/*{\Mop returns one edge with source |v|. It's the starting point for
  the circular iteration over the edges with source |v|.
  \precond |!is_isolated(v)|.}*/
{ return v->out_sedge_; }

SHalfedge_const_handle last_out_edge(SVertex_const_handle v) const
/*{\Mop returns one edge with source |v|. \precond |!is_isolated(v)|.}*/
{ return cap(v->out_sedge_); }

SHalfedge_const_handle cyclic_adj_succ(SHalfedge_const_handle e) const
/*{\Mop returns the edge after |e| in the cyclic ordered adjacency list of
  |source(e)|.}*/
{ return e->sprev_->twin_; }

SHalfedge_const_handle cyclic_adj_pred(SHalfedge_const_handle e) const
/*{\Mop returns the edge before |e| in the cyclic ordered adjacency list of
  |source(e)|.}*/
{ return e->twin_->snext_; }


SHalfedge_const_handle next(SHalfedge_const_handle e) const
/*{\Mop returns the next edge in the face cycle containing |e|.}*/
{ return e->snext_; }

SHalfedge_const_handle previous(SHalfedge_const_handle e) const
/*{\Mop returns the previous edge in the face cycle containing |e|.}*/
{ return e->sprev_; }

SFace_const_handle face(SHalfedge_const_handle e) const
/*{\Mop returns the face incident to |e|.}*/
{ return e->incident_sface_; }

SFace_const_handle face(SHalfloop_const_handle l) const
/*{\Mop returns the face incident to |l|.}*/
{ return l->incident_sface_; }

SFace_const_handle face(SVertex_const_handle v) const
/*{\Mop returns the face incident to |v|.
   \precond |is_isolated(v)|.}*/
{ return v->incident_sface_; }

/*{\Mtext \headerline{Iteration} \setopdims{3.3cm}{0cm}}*/
  
SVertex_const_iterator svertices_begin() const
{ return psm_->svertices_begin(); }
SVertex_const_iterator svertices_end() const
{ return psm_->svertices_end(); }
SHalfedge_const_iterator shalfedges_begin() const
{ return psm_->shalfedges_begin(); }
SHalfedge_const_iterator shalfedges_end() const
{ return psm_->shalfedges_end(); }
SFace_const_iterator sfaces_begin() const
{ return psm_->sfaces_begin(); }
SFace_const_iterator sfaces_end() const
{ return psm_->sfaces_end(); }
SHalfloop_const_iterator shalfloops_begin() const
{ return psm_->shalfloops_begin(); }
SHalfloop_const_iterator shalfloops_end() const
{ return psm_->shalfloops_end(); }

SHalfloop_const_handle shalfloop() const
/*{\Mop returns access to the loop.}*/
{ return psm_->shalfloop(); }

bool has_loop() const
/*{\Mop returns true iff there is a loop.}*/
{ return shalfloop() != sncp()->shalfloops_end(); }

SHalfedge_around_svertex_const_circulator 
  out_edges(SVertex_const_handle v) const
/*{\Mop returns a circulator for the cyclic adjacency list of |v|.
\precond the adjacency list is not empty.}*/
{ return SHalfedge_around_svertex_const_circulator(first_out_edge(v)); }

SFace_cycle_const_iterator sface_cycles_begin(SFace_const_handle f) const
/*{\Mop returns an iterator for all bounding face cycles of |f|.
The iterator is is convertable to |SVertex_const_handle|, 
|SHalfloop_const_handle|, or |SHalfedge_const_handle|.}*/
{ return f->sface_cycles_begin(); }

SFace_cycle_const_iterator sface_cycles_end(SFace_const_handle f) const
/*{\Mop returns the past the end iterator of |f|.}*/
{ return f->sface_cycles_end(); }

/*{\Mtext \headerline{Statistics and Integrity}}*/

Size_type number_of_svertices() const 
/*{\Mop returns the number of vertices.}*/
{ Size_type n(0);
  SVertex_const_iterator vit;
  CGAL_nef3_forall_svertices_of(vit,psm_) ++n;
  return n; }

Size_type number_of_shalfedges() const 
/*{\Mop returns the number of halfedges.}*/
{ Size_type n(0);
  SHalfedge_const_iterator eit;
  CGAL_nef3_forall_shalfedges_of(eit,psm_) ++n;
  return n; }

Size_type number_of_sedges() const 
/*{\Mop returns the number of edges.}*/
{ return number_of_shalfedges()/2; }

Size_type number_of_shalfloops() const 
/*{\Mop returns the number of halfloops.}*/
{ return ( has_loop() ? 2 : 0); }

Size_type number_of_sloops() const 
/*{\Mop returns the number of loops.}*/
{ return number_of_shalfloops()/2; }

Size_type number_of_sfaces() const    
/*{\Mop returns the number of faces.}*/
{ Size_type n(0);
  SFace_const_iterator fit;
  CGAL_nef3_forall_sfaces_of(fit,psm_) ++n;
  return n; }

Size_type number_of_sface_cycles() const;
/*{\Mop returns the number of non-trivial sface cycles.}*/

Size_type number_of_connected_components() const;
/*{\Mop calculates the number of connected components of |P|.}*/

void print_statistics(std::stringstream& os = std::cout) const
/*{\Mop print the statistics of |P|: the number of vertices, edges, 
and faces.}*/
{
  os << "Sphere Map - Statistics\n";
  os << "|V| = " << number_of_svertices() << std::endl;
  os << "|E| = " << number_of_shalfedges() << std::endl;
  os << "|L| = " << number_of_shalfloops() << std::endl;
  os << "|F| = " << number_of_sfaces() << std::endl;
  os << "|Fcs| = " << number_of_sface_cycles() << std::endl << std::endl;
}
 
void check_integrity_and_topological_planarity(bool faces=true) const;
/*{\Mop checks the link structure and the genus of |P|.}*/

SHalfedge_const_handle cas(SHalfedge_const_handle e) const 
{ return cyclic_adj_succ(e); } 

SHalfedge_const_handle cap(SHalfedge_const_handle e) const
{ return cyclic_adj_pred(e); }

template <typename H>
bool is_boundary_object(H h) const
{ return psm_->is_sm_boundary_object(h); }

/*{\Mtext \headerline{Associated Information}\restoreopdims}*/

const Sphere_point& point(SVertex_const_handle v) const
/*{\Mop returns the embedding of |v|.}*/
{ return v->tmp_point(); }

const Sphere_circle& circle(SHalfedge_const_handle e) const
/*{\Mop returns the circle supporting |e|.}*/
{ return e->tmp_circle(); }

const Sphere_circle& circle(SHalfloop_const_handle l) const
/*{\Mop returns the circle supporting |l|.}*/
{ return l->tmp_circle(); }

Mark mark(SVertex_const_handle v) const
/*{\Mop returns the mark of |v|.}*/
{ return v->mark_; }

Mark mark(SHalfedge_const_handle e) const
/*{\Mop returns the mark of |e|.}*/
{ return ( &*e < &*twin(e) ) ? e->mark_ : twin(e)->mark_; }

Mark mark(SHalfloop_const_handle l) const
/*{\Mop returns the mark of |l|.}*/
{ return ( &*l < &*twin(l) ) ? l->mark_ : twin(l)->mark_; }

Mark mark(SFace_const_handle f) const
/*{\Mop returns the mark of |f|.}*/
{ return f->mark_; }

Mark mark_of_halfsphere(int i) const
{ CGAL_nef3_assertion(i);
  if (i<0) return psm_->m_neg_; 
  return psm_->m_pos_; }

/*{\Mtext \headerline{Iteration}}*/
/*{\Mtext The list of all objects can be accessed via iterator ranges.
For comfortable iteration we also provide iterations macros. 
The iterator range access operations are of the following kind:\\
|SVertex_iterator   svertices_begin()/svertices_end()|\\
|SHalfedge_iterator shalfedges_begin()/shalfedges_end()|\\
|SHalfloop_iterator shalfloops_begin()/shalfloops_end()|\\
|SFace_iterator     sfaces_begin()/sfaces_end()|

The macros are then |CGAL_nef3_forall_svertices(v,M)|,
|CGAL_nef3_forall_shalfedges(e,M)|, |CGAL_nef3_forall_sedges(e,M)|,
|CGAL_nef3_forall_sfaces(f,M)|, |CGAL_nef3_forall_sface_cycles_of(fc,F)| 
where |M|is a sphere map and |F| is a face.}*/

}; // SNC_SM_const_decorator


template <typename SM_>
void SNC_SM_const_decorator<SM_>::
check_integrity_and_topological_planarity(bool faces) const
{
  TRACEN("check_integrity_and_topological_planarity:");
  using CGAL::Object_index;
  Object_index<SVertex_const_iterator>   
    VI(svertices_begin(),svertices_end(),'v');
  Object_index<SHalfedge_const_iterator> 
    EI(shalfedges_begin(),shalfedges_end(),'e');
  Object_index<SFace_const_iterator> 
    FI(sfaces_begin(),sfaces_end(),'f');
  typedef SHalfedge_around_svertex_const_circulator hvc_circulator;
  typedef SHalfedge_around_sface_const_circulator   hfc_circulator;
  SVertex_const_handle v;
  int iso_vert_num=0;
  /* check the source links of out edges and count isolated vertices */
  CGAL_nef3_forall_svertices(v,*this) {
    if ( is_isolated(v) ) {
      if ( faces )
        CGAL_nef3_assertion_msg(face(v) != SFace_const_handle(), 
				VI(v).c_str());
      ++iso_vert_num;
    } else {
      CGAL_nef3_assertion_msg(first_out_edge(v) != SHalfedge_const_handle(),
      VI(v).c_str());
      TRACEN(point(v)<<" "<<EI(first_out_edge(v)));
      CGAL_nef3_assertion_msg(source(first_out_edge(v)) == v ,
			      VI(v).c_str());
    }
  }

  /* check the bidirected links and the face pointer init */
  SHalfedge_const_iterator e;
  CGAL_nef3_forall_shalfedges(e,*this) {
    CGAL_nef3_assertion( twin(twin(e)) == e );
    CGAL_nef3_assertion( source(e) != SVertex_const_handle() );
    CGAL_nef3_assertion( next(e) != SHalfedge_const_handle() );
    CGAL_nef3_assertion( previous(next(e)) == e );
    CGAL_nef3_assertion( target(e) == source(next(e)) );
    CGAL_nef3_assertion( previous(e) != SHalfedge_const_handle() );
    CGAL_nef3_assertion( next(previous(e)) == e );
    CGAL_nef3_assertion( target(previous(e)) == source(e) );
    if ( !faces ) continue;
    CGAL_nef3_assertion( face(e) != SFace_const_handle() );
    CGAL_nef3_assertion( face(next(e)) == face(e) );
    CGAL_nef3_assertion( face(previous(e)) == face(e) );
  }

  int fc_num(0),iv_num(0);
  SFace_const_iterator f;
  SFace_cycle_const_iterator fci;
  CGAL_nef3_forall_sfaces(f,*this) {
    CGAL_nef3_forall_sface_cycles_of(fci,f) {
      if ( fci.is_shalfedge() ) {
        CGAL_nef3_assertion( face(SHalfedge_const_handle(fci)) == f ); 
	++fc_num;
      } else if ( fci.is_svertex() ) {
        CGAL_nef3_assertion( face(SVertex_const_handle(fci)) == f ); 
	++iv_num;
      } else if( fci.is_shalfloop() ) {
        CGAL_nef3_assertion( face(SHalfloop_const_handle(fci)) == f );
	++fc_num;
      }
      else CGAL_nef3_assertion(0);
    }
  }

  int v_num = number_of_svertices() - iso_vert_num + 
    number_of_shalfloops();
  int e_num = number_of_sedges() + 
    number_of_shalfloops();
  int c_num = number_of_connected_components() - iso_vert_num + 
    number_of_sloops();
  int f_num = number_of_sface_cycles() - c_num + 1;
  TRACEV(fc_num);TRACEV(iv_num);TRACEV(iso_vert_num);
  TRACEV(v_num);TRACEV(e_num);TRACEV(c_num);TRACEV(f_num);
  /* this means all face cycles and all isolated vertices are 
     indeed referenced from a face */
  /* every isolated vertex increases the component count
       one face cycle per component is redundent except one
       finally check the Euler formula: */
  CGAL_nef3_assertion( v_num - e_num + f_num == 1 + c_num );
}

template <typename SM_>
typename SNC_SM_const_decorator<SM_>::Size_type
SNC_SM_const_decorator<SM_>::
number_of_sface_cycles() const
{
  unsigned int fc_num=0;
  CGAL::Unique_hash_map<SHalfedge_const_handle,bool> visited; 
    // init with bool() == false
  SHalfedge_const_iterator e;
  CGAL_nef3_forall_shalfedges(e,*this) {
    if ( visited[e] ) continue;
    SHalfedge_around_sface_const_circulator hfc(e), hend(hfc);
    CGAL_For_all(hfc,hend) visited[hfc]=true;
    ++fc_num;
  }
  if ( has_loop() ) fc_num += 2;
  return fc_num;
}

template <typename SM_>
typename SNC_SM_const_decorator<SM_>::Size_type
SNC_SM_const_decorator<SM_>::
number_of_connected_components() const
{
  int comp_num=0;
  CGAL::Unique_hash_map<SVertex_const_iterator,bool> visited(false); 
  SVertex_const_iterator v;
  CGAL_nef3_forall_svertices(v,*this) {
    if (visited[v]) continue;
    std::list<SVertex_const_iterator> L;
    L.push_back(v); visited[v]=true; 
    /* we keep the invariant that all nodes which have been stacked
       are marked visited */
    while (!L.empty()) {
      SVertex_const_iterator vc = L.front(); L.pop_front();
      if ( is_isolated(vc) ) continue;
      SHalfedge_around_svertex_const_circulator 
	havc(first_out_edge(vc)), hend(havc);
      CGAL_For_all(havc,hend) {
        if (!visited[target(havc)]) {
          L.push_back(target(havc)); visited[target(havc)]=true; 
        }
      }
    }
    ++comp_num;
  }
  return comp_num;   
}

CGAL_END_NAMESPACE
#endif // CGAL_SNC_SM_CONST_DECORATOR_H

