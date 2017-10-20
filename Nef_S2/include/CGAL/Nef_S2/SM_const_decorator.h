// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de>
//                 Miguel Granados    <granados@mpi-sb.mpg.de>
//                 Susan Hert         <hert@mpi-sb.mpg.de>
//                 Lutz Kettner       <kettner@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_SM_CONST_DECORATOR_H 
#define CGAL_SM_CONST_DECORATOR_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/basic.h>
#include <CGAL/circulator.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_2/Object_index.h>
#include <CGAL/Nef_S2/SM_iteration.h>
#include <CGAL/Nef_S2/SM_decorator_traits.h>
#include <string>
#include <list>
#include <sstream>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 67
#include <CGAL/Nef_2/debug.h>
#ifndef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <boost/any.hpp>
#endif

namespace CGAL {

template <typename Map_>
class SM_const_decorator { 

  typedef SM_const_decorator<Map_> Self;
public:

  typedef Map_                       Map;
  typedef const Map_                 Sphere_map;
  typedef SM_decorator_const_traits<Map>  Decorator_traits;

/*{\Mdefinition ...}*/

/*{\Mtypes 5}*/

typedef typename Map::Sphere_kernel Sphere_kernel;
/*{\Mtypemember spherical geometry.}*/

typedef typename Map::Sphere_point   Sphere_point;
/*{\Mtypemember embedding vertices.}*/

typedef typename Map::Sphere_segment Sphere_segment;
/*{\Mtypemember embedding edges.}*/

typedef typename Map::Sphere_circle  Sphere_circle;
/*{\Mtypemember embedding loops.}*/

typedef typename Map::Sphere_direction Sphere_direction;
/*{\Mtypemember embedding directions.}*/

typedef typename Map::Mark   Mark;
/*{\Mtypemember attributes of objects (vertices, edges, faces).}*/

typedef size_t Size_type;
/*{\Mtypemember size type.}*/

#ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
typedef void*  GenPtr;
#else
typedef boost::any GenPtr;
#endif


// typedef typename Map::Constructor_const_parameter Constructor_parameter;
typedef typename Map::SVertex_const_handle SVertex_const_handle;
typedef typename Map::SVertex_const_iterator SVertex_const_iterator;
typedef typename Map::SHalfedge_const_handle SHalfedge_const_handle; 
typedef typename Map::SHalfedge_const_iterator SHalfedge_const_iterator;
typedef typename Map::SHalfloop_const_handle SHalfloop_const_handle;
typedef typename Map::SHalfloop_const_iterator SHalfloop_const_iterator;
typedef typename Map::SFace_const_handle SFace_const_handle;
typedef typename Map::SFace_const_iterator SFace_const_iterator;

/*{\Mtext Local types are handles, iterators and circulators of the
following kind: |SVertex_handle|, |SVertex_iterator|, |SHalfedge_handle|,
|SHalfedge_iterator|, |SHalfloop_handle|, |SHalfloop_iterator|,
|SFace_handle|, |SFace_iterator|.  Additionally the following
circulators are defined.}*/

typedef typename Map::SHalfedge_around_svertex_const_circulator
        SHalfedge_around_svertex_const_circulator;
/*{\Mtypemember circulating the adjacency list of an vertex |v|.}*/

typedef typename Map::SHalfedge_around_sface_const_circulator
        SHalfedge_around_sface_const_circulator;
/*{\Mtypemember circulating the face cycle of an face |f|.}*/

typedef typename Map::SFace_cycle_const_iterator 
        SFace_cycle_const_iterator;
/*{\Mtypemember iterating all sface cycles of an sface |f|.
The iterator has method |bool is_svertex()|, |bool is_shalfedge()|,
|bool is_shalfloop()|, and can be converted to the corresponding
handles |SVertex_const_handle|, |SHalfedge_const_handle|, or 
|SHalfloop_const_handle|.}*/

protected: 
  const Map* psm_;

  void set_sm(const Map* W) {
    psm_ = W;
  }

public:

/*{\Mcreation 3}*/
SM_const_decorator() : psm_(0) {}
SM_const_decorator(const Self& D) : psm_(D.psm_) {}
Self& operator=(const Self& D) { psm_=D.psm_; return *this; }

SM_const_decorator(const Map* M) : psm_(M) {}
/*{\Mcreate constructs a plane map decorator exploring |M|.}*/

/*{\Moperations 4 4}*/

const Map* sphere_map() const { return psm_; }

bool is_isolated(SVertex_const_handle v) const
/*{\Mop returns |true| when |v| is linked to the interior of a face.}*/
{ return (SHalfedge_const_handle(v->out_sedge()) == SHalfedge_const_handle()); }

SHalfedge_const_handle first_out_edge(SVertex_const_handle v) const
/*{\Mop returns one edge with source |v|. It's the starting point for
  the circular iteration over the edges with source |v|.
  \precond |!is_isolated(v)|.}*/
{ return v->out_sedge(); }

SHalfedge_const_handle last_out_edge(SVertex_const_handle v) const
/*{\Mop returns one edge with source |v|. \precond |!is_isolated(v)|.}*/
{ return cap(v->out_sedge()); }

SHalfedge_const_handle cyclic_adj_succ(SHalfedge_const_handle e) const
/*{\Mop returns the edge after |e| in the cyclic ordered adjacency list of
  |e->source()|.}*/
{ return e->sprev()->twin(); }

SHalfedge_const_handle cyclic_adj_pred(SHalfedge_const_handle e) const
/*{\Mop returns the edge before |e| in the cyclic ordered adjacency list of
  |e->source()|.}*/
{ return e->twin()->snext(); }


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

bool has_shalfloop() const
/*{\Mop returns true iff there is a loop.}*/
{ return psm_->has_shalfloop(); }

SHalfedge_around_svertex_const_circulator 
  out_edges(SVertex_const_handle v) const
/*{\Mop returns a circulator for the cyclic adjacency list of |v|.
\precond the adjacency list is not empty.}*/
{ return SHalfedge_around_svertex_const_circulator(first_out_edge(v)); }

SFace_cycle_const_iterator sface_cycles_begin(SFace_const_handle f) const
/*{\Mop returns an iterator for all bounding face cycles of |f|.
The iterator is is convertable to |SVertex_const_handle|, 
|SHalfloop_const_handle|, or |SHalfedge_const_handle|.}*/
{ return f->boundary_entry_objects_.begin(); }

SFace_cycle_const_iterator sface_cycles_end(SFace_const_handle f) const
/*{\Mop returns the past the end iterator of |f|.}*/
{ return f->boundary_entry_objects_.end(); }

/*{\Mtext \headerline{Statistics and Integrity}}*/

Size_type number_of_svertices() const 
/*{\Mop returns the number of vertices.}*/
{ return psm_->number_of_svertices(); }

Size_type number_of_shalfedges() const 
/*{\Mop returns the number of halfedges.}*/
{ return psm_->number_of_shalfedges(); }

Size_type number_of_sedges() const 
/*{\Mop returns the number of edges.}*/
{ return number_of_shalfedges()/2; }

Size_type number_of_shalfloops() const 
/*{\Mop returns the number of halfloops.}*/
{ return psm_->number_of_shalfloops(); }

Size_type number_of_sloops() const 
/*{\Mop returns the number of loops.}*/
{ return psm_->number_of_shalfloops()/2; }

Size_type number_of_sfaces() const    
/*{\Mop returns the number of faces.}*/
{ return psm_->number_of_sfaces(); }

Size_type number_of_sface_cycles() const;
/*{\Mop returns the number of face cycles.}*/

Size_type number_of_connected_components() const;
/*{\Mop calculates the number of connected components of |P|.}*/

void print_statistics(std::ostream& os = std::cout) const
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
bool is_sm_boundary_object(H h) const
{ return psm_->is_sm_boundary_object(h); }

/*{\Mtext \headerline{Associated Information}\restoreopdims}*/

/*{\Mtext \headerline{Iteration}}*/
/*{\Mtext The list of all objects can be accessed via iterator ranges.
For comfortable iteration we also provide iterations macros. 
The iterator range access operations are of the following kind:\\
|SVertex_iterator   svertices_begin()/svertices_end()|\\
|SHalfedge_iterator shalfedges_begin()/shalfedges_end()|\\
|SHalfloop_iterator shalfloops_begin()/shalfloops_end()|\\
|SFace_iterator     sfaces_begin()/sfaces_end()|

The macros are then |CGAL_forall_svertices(v,M)|,
|CGAL_forall_shalfedges(e,M)|, |CGAL_forall_sedges(e,M)|,
|CGAL_forall_sfaces(f,M)|, |CGAL_forall_sface_cycles_of(fc,F)|
where |M| is a sphere map and |F| is a sface.}*/

}; // SM_const_decorator



template <typename SM_>
void SM_const_decorator<SM_>::
check_integrity_and_topological_planarity(bool faces) const
{
  CGAL_NEF_TRACEN("check_integrity_and_topological_planarity:");
  using CGAL::Object_index;
  Object_index<SVertex_const_iterator>   
    VI(svertices_begin(),svertices_end(),'v');
  Object_index<SHalfedge_const_iterator> 
    EI(shalfedges_begin(),shalfedges_end(),'e');
  Object_index<SFace_const_iterator> 
    FI(sfaces_begin(),sfaces_end(),'f');
  SVertex_const_handle v;
  int iso_vert_num=0;
  /* check the source links of out edges and count isolated vertices */
  CGAL_forall_svertices(v,*this) {
    if ( is_isolated(v) ) {
      if ( faces )
        CGAL_assertion_msg(v->incident_sface() != SFace_const_handle(), VI(v).c_str());
      ++iso_vert_num;
    } else {
      CGAL_assertion_msg(first_out_edge(v) != SHalfedge_const_handle(),
      VI(v).c_str());
      CGAL_NEF_TRACEN(v->point()<<" "<<EI(first_out_edge(v)));
      CGAL_assertion_msg(first_out_edge(v)->source() == v,
			 VI(v).c_str());
    }
  }

  /* check the bidirected links and the face pointer init */
  SHalfedge_const_iterator e;
  CGAL_forall_shalfedges(e,*this) {
    CGAL_assertion( e->twin()->twin() == e );
    CGAL_assertion( e->source() != SVertex_const_handle() );
    CGAL_assertion( e->snext() != SHalfedge_const_handle() );
    CGAL_assertion( e->snext()->sprev() == e );
    CGAL_assertion( e->target() == e->snext()->source() );
    CGAL_assertion( e->sprev() != SHalfedge_const_handle() );
    CGAL_assertion( e->sprev()->snext() == e );
    CGAL_assertion( e->sprev()->twin()->source() == e->source() );
    if ( !faces ) continue;
    CGAL_assertion( e->incident_sface() != SFace_const_handle() );
    CGAL_assertion( e->snext()->incident_sface() == e->incident_sface() );
    CGAL_assertion( e->sprev()->incident_sface() == e->incident_sface() );
  }

  int fc_num(0),iv_num(0);
  SFace_const_iterator f;
  SFace_cycle_const_iterator fci;
  CGAL_forall_sfaces(f,*this) {
    CGAL_forall_sface_cycles_of(fci,f) {
      if ( fci.is_shalfedge() ) {
        CGAL_assertion( SHalfedge_const_handle(fci)->incident_sface() == f ); 
	++fc_num;
      } else if ( fci.is_svertex() ) {
        CGAL_assertion( SVertex_const_handle(fci)->incident_sface() == f ); 
	++iv_num;
      } else if( fci.is_shalfloop() ) {
        CGAL_assertion( SHalfloop_const_handle(fci)->incident_sface() == f );
	++fc_num;
      } else CGAL_error_msg("damn generic handle.");
    }
  }

  CGAL_assertion_code(std::size_t v_num = number_of_svertices() - 
		      iso_vert_num + 
		      number_of_shalfloops());
  CGAL_assertion_code(std::size_t e_num = number_of_sedges() + 
		      number_of_shalfloops());
  CGAL_assertion_code(std::size_t c_num = number_of_connected_components() - 
		      iso_vert_num 
		      + number_of_sloops());
  CGAL_assertion_code(std::size_t f_num = number_of_sface_cycles() - c_num + 1);
  CGAL_assertion_code(CGAL_NEF_TRACEV(fc_num));
  CGAL_assertion_code(CGAL_NEF_TRACEV(iv_num));
  CGAL_assertion_code(CGAL_NEF_TRACEV(iso_vert_num));
  CGAL_assertion_code(CGAL_NEF_TRACEV(v_num));
  CGAL_assertion_code(CGAL_NEF_TRACEV(e_num));
  CGAL_assertion_code(CGAL_NEF_TRACEV(c_num));
  CGAL_assertion_code(CGAL_NEF_TRACEV(f_num));
  /* this means all face cycles and all isolated vertices are 
     indeed referenced from a face */
  /* every isolated vertex increases the component count
       one face cycle per component is redundent except one
       finally check the Euler formula: */
  CGAL_assertion( v_num - e_num + f_num == 1 + c_num );
}

template <typename SM_>
typename SM_const_decorator<SM_>::Size_type
SM_const_decorator<SM_>::
number_of_sface_cycles() const
{
  unsigned int fc_num=0;
  CGAL::Unique_hash_map<SHalfedge_const_handle,bool> visited; 
  SHalfedge_const_iterator e;
  CGAL_forall_shalfedges(e,*this) {
    if (visited[e]) continue;
    SHalfedge_around_sface_const_circulator hfc(e), hend(hfc);
    CGAL_For_all(hfc,hend) visited[hfc]=true;
    ++fc_num;
  }
  if ( has_shalfloop() ) fc_num += 2;
  return fc_num;
}

template <typename SM_>
typename SM_const_decorator<SM_>::Size_type
SM_const_decorator<SM_>::
number_of_connected_components() const
{
  int comp_num=0;
  CGAL::Unique_hash_map<SVertex_const_iterator,bool> visited(false); 
  SVertex_const_iterator v;
  CGAL_forall_svertices(v,*this) {
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
        if (!visited[havc->target()]) {
          L.push_back(havc->target()); visited[havc->target()]=true; 
        }
      }
    }
    ++comp_num;
  }
  return comp_num;   
}

} //namespace CGAL
#endif // CGAL_SM_CONST_DECORATOR_H
