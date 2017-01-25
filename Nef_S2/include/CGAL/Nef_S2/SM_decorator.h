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
// 
//
// Author(s)     : Michael Seel       <seel@mpi-sb.mpg.de>
//                 Miguel Granados    <granados@mpi-sb.mpg.de>
//                 Susan Hert         <hert@mpi-sb.mpg.de>
//                 Lutz Kettner       <kettner@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_SM_DECORATOR_H 
#define CGAL_SM_DECORATOR_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/basic.h>
#include <CGAL/Nef_S2/SM_const_decorator.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG  23
#include <CGAL/Nef_2/debug.h>
#include <CGAL/Nef_S2/SM_decorator_traits.h>
#include <CGAL/Nef_S2/Sphere_map.h>
#include <CGAL/Nef_S2/Normalizing.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/IO/Verbose_ostream.h>
#ifndef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <boost/any.hpp>
#endif

namespace CGAL {

/*{\Moptions print_title=yes }*/ 
/*{\Moptions outfile=SM_decorator.man }*/
/*{\Manpage {SM_decorator}{Sphere_map}
{Topological sphere map decorator}{D}}*/

template <typename Map_>
class SM_decorator
{ 
public:
typedef Map_              Map;
typedef Map_              Sphere_map;
typedef SM_decorator<Map> Self;
typedef SM_decorator_traits<Map>  Decorator_traits;
/*{\Mdefinition ...}*/

/*{\Mtypes 5}*/

typedef CGAL::SM_const_decorator<Map>   SM_const_decorator;

typedef typename Map::Sphere_kernel    Sphere_kernel;
/*{\Mtypemember spherical geometry.}*/

typedef typename Map::Sphere_point     Sphere_point;
/*{\Mtypemember embedding vertices.}*/

typedef typename Map::Sphere_segment   Sphere_segment;
/*{\Mtypemember embedding edges.}*/

typedef typename Map::Sphere_circle    Sphere_circle;
/*{\Mtypemember embedding loops.}*/

typedef typename Map::Sphere_direction Sphere_direction;
/*{\Mtypemember embedding directions.}*/

typedef typename Map::Mark   Mark;
/*{\Mtypemember attributes of objects (vertices, edges, faces).}*/

typedef size_t Size_type;
/*{\Mtypemember size type.}*/

enum { BEFORE = -1, AFTER = 1 };
/*{\Menum insertion order labels.}*/

typedef typename Sphere_kernel::Aff_transformation_3 Aff_transformation_3;

#ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
typedef void*  GenPtr;
#else
typedef boost::any GenPtr;
#endif
typedef typename Map::SVertex                   SVertex;
typedef typename Map::SVertex_handle            SVertex_handle;
typedef typename Map::SVertex_iterator          SVertex_iterator;
typedef typename Map::SVertex_const_handle      SVertex_const_handle;
typedef typename Map::SVertex_const_iterator    SVertex_const_iterator;
typedef typename Map::SHalfedge                 SHalfedge;
typedef typename Map::SHalfedge_handle          SHalfedge_handle;
typedef typename Map::SHalfedge_iterator        SHalfedge_iterator;
typedef typename Map::SHalfedge_const_handle    SHalfedge_const_handle;
typedef typename Map::SHalfedge_const_iterator  SHalfedge_const_iterator;
typedef typename Map::SFace                     SFace;
typedef typename Map::SFace_handle              SFace_handle;
typedef typename Map::SFace_iterator            SFace_iterator;
typedef typename Map::SFace_const_handle        SFace_const_handle;
typedef typename Map::SFace_const_iterator      SFace_const_iterator;
typedef typename Map::SHalfloop                 SHalfloop;
typedef typename Map::SHalfloop_handle          SHalfloop_handle;
typedef typename Map::SHalfloop_iterator        SHalfloop_iterator;
typedef typename Map::SHalfloop_const_handle    SHalfloop_const_handle;
typedef typename Map::SHalfloop_const_iterator  SHalfloop_const_iterator;
typedef typename Map::Object_handle             Object_handle;
typedef typename Map::SHalfedge_around_svertex_circulator
                      SHalfedge_around_svertex_circulator;
typedef typename Map::SHalfedge_around_sface_circulator
                      SHalfedge_around_sface_circulator;
typedef typename Map::SFace_cycle_iterator      SFace_cycle_iterator;

/*{\Mtypemember iterating all face cycles of an face |f|.
The iterator has method |bool is_svertex()|, |bool is_shalfedge()|,
|bool is_shalfloop()|, and can be converted to the corresponding
handles |SVertex_handle|, |SHalfedge_handle|, or 
|SHalfloop_handle|.}*/

protected: 
  // don't change this into a shared_ptr even if it seems sensible.
  // minkowski_sum_3 already has a fix in place that deletes the
  // object psm_ points to.
  Map* psm_;

public:
/*{\Mcreation 3}*/
SM_decorator() : psm_(0) {}
SM_decorator(const Self& D) : psm_(D.psm_) {}
Self& operator=(const Self& D) { psm_=D.psm_; return *this; }

SM_decorator(Map* M) : psm_(M) {}
 
/*{\Moperations 4 4}*/

Map* sphere_map() const { return psm_; }

Map* map() { return psm_; }
const Map* map() const { return psm_; }

bool is_isolated(SVertex_const_handle v) const
{ return (v->out_sedge() == SHalfedge_handle()); }

bool is_isolated(SVertex_handle v) const
/*{\Mop returns |true| when |v| is linked to the interior of a face.}*/
{ return (v->out_sedge() == SHalfedge_handle()); }

SHalfedge_const_handle first_out_edge(SVertex_const_handle v) const
{ return v->out_sedge(); }
SHalfedge_const_handle last_out_edge(SVertex_const_handle v) const
{ return cap(v->out_sedge()); }

SHalfedge_handle first_out_edge(SVertex_handle v) const
/*{\Mop returns one edge with source |v|. It's the starting point for
  the circular iteration over the edges with source |v|.
  \precond |!is_isolated(v)|.}*/
{ return v->out_sedge(); }

SHalfedge_handle last_out_edge(SVertex_handle v) const
/*{\Mop returns one edge with source |v|. \precond |!is_isolated(v)|.}*/
{ return cap(v->out_sedge()); }

SHalfedge_const_handle cyclic_adj_succ(SHalfedge_const_handle e) const
{ return e->sprev()->twin(); }
SHalfedge_const_handle cyclic_adj_pred(SHalfedge_const_handle e) const
{ return e->twin()->snext(); }

SHalfedge_handle cyclic_adj_succ(SHalfedge_handle e) const
/*{\Mop returns the edge after |e| in the cyclic ordered adjacency list of
  |e->source()|.}*/
{ return e->sprev()->twin(); }

SHalfedge_handle cyclic_adj_pred(SHalfedge_handle e) const
/*{\Mop returns the edge before |e| in the cyclic ordered adjacency list of
  |e->source()|.}*/
{ return e->twin()->snext(); }

bool has_shalfloop() const
/*{\Mop returns true iff there is a loop.}*/
{ return psm_->has_shalfloop(); }

/*{\Mtext \headerline{Iteration} \setopdims{3.3cm}{0cm}}*/
  
SVertex_iterator svertices_begin() const
{ return psm_->svertices_begin(); }
SVertex_iterator svertices_end() const
{ return psm_->svertices_end(); }
SHalfedge_iterator shalfedges_begin() const
{ return psm_->shalfedges_begin(); }
SHalfedge_iterator shalfedges_end() const
{ return psm_->shalfedges_end(); }
SFace_iterator sfaces_begin() const
{ return psm_->sfaces_begin(); }
SFace_iterator sfaces_end() const
{ return psm_->sfaces_end(); }
SHalfloop_iterator shalfloops_begin() const
{ return psm_->shalfloops_begin(); }
SHalfloop_iterator shalfloops_end() const
{ return psm_->shalfloops_end(); }

SHalfloop_handle& shalfloop()
{ return psm_->shalfloop(); }
SHalfloop_const_handle shalfloop() const
{ return psm_->shalfloop(); }

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

SFace_cycle_iterator sface_cycles_begin(SFace_handle f) const
/*{\Mop returns an iterator for all bounding face cycles of |f|.
The iterator is is convertable to |SVertex_handle|, 
|SHalfloop_handle|, or |SHalfedge_handle|.}*/
{ return f->boundary_entry_objects().begin(); }

SFace_cycle_iterator sface_cycles_end(SFace_handle f) const
/*{\Mop returns the past the end iterator of |f|.}*/
{ return f->boundary_entry_objects_.end(); }

SHalfedge_around_svertex_circulator 
  out_edges(SVertex_handle v) const
/*{\Mop returns a circulator for the cyclic adjacency list of |v|.
\precond the adjacency list is not empty.}*/
{ return SHalfedge_around_svertex_circulator(first_out_edge(v)); }

/*{\Mtext \headerline{Update Operations}}*/

void clear() const
/*{\Mop reintializes |P| to the empty map.}*/
{ psm_->clear(); }

bool is_closed_at_source(SHalfedge_handle e) const
/*{\Mop returns |true| when |e->sprev() == e->twin()|.}*/
{ return e->sprev() == e->twin(); }

bool is_closed_at_target(SHalfedge_handle e) const
/*{\Mop returns |true| when |e->snext() == e->twin()|.}*/
{ return e->snext() == e->twin(); }

SHalfedge_handle cas(SHalfedge_handle e) const 
{ return cyclic_adj_succ(e); } 
SHalfedge_handle cap(SHalfedge_handle e) const
{ return cyclic_adj_pred(e); }

template <typename H>
void make_twins(H h1, H h2) const
{ h1->twin() = h2; h2->twin() = h1; }

SVertex_handle new_svertex(const Sphere_point& p = Sphere_point())
/*{\Mop creates a new vertex.}*/
{ return map()->new_svertex(p); }

SHalfedge_handle new_shalfedge_pair() {
/*{\Xop creates a new edge pair. No connectivity is provided.}*/
  return map()->new_shalfedge_pair();
}

SHalfloop_handle new_shalfloop_pair()
/*{\Mop creates a new loop pair.
\precond No sloop pair exists in the local graph.}*/ 
{ CGAL_assertion(!has_shalfloop());
  return map()->new_shalfloop_pair(); 
}

SFace_handle new_sface()
/*{\Mop creates a new face.}*/
{ return map()->new_sface(); }

void delete_vertex_only(SVertex_handle v)
/*{\Mop deletes |v| without any connectivity update.}*/
{ map()->delete_svertex(v); }

void delete_edge_pair_only(SHalfedge_handle e)
/*{\Mop deletes |e| and its twin without any connectivity update.}*/
{ map()->delete_shalfedge_pair(e); }

void delete_halfedge_only(SHalfedge_handle e)
/*{\Mop deletes |e| without its twin and without any connectivity update.}*/
{ map()->delete_shalfedge(e); }

void delete_face_only(SFace_handle f) 
/*{\Mop deletes |f| without any connectivity update.}*/
{ map()->delete_sface(f); }

void delete_loop_only()
/*{\Mop deletes the loop and its twin without any connectivity update.}*/ 
{ map()->delete_shalfloop_pair(); }

template <typename H>
bool is_sm_boundary_object(H h) const
{ return map()->is_sm_boundary_object(h); }

template <typename H>
void store_sm_boundary_object(H h, SFace_handle f) {
  CGAL_assertion(!map()->is_sm_boundary_object(h));
  f->boundary_entry_objects().push_back(make_object(h));
  map()->store_sm_boundary_item(h, --(f->sface_cycles_end()));
}

template <typename H>
void undo_sm_boundary_object(H h, SFace_handle f)
{ CGAL_assertion(map()->is_sm_boundary_object(h));
  SFace_cycle_iterator it = map()->sm_boundary_item(h);
  map()->undef_sm_boundary_item(h);
  f->boundary_entry_objects().erase(it);
}

void link_as_face_cycle(SHalfedge_handle e, SFace_handle f)
/*{\Mop creates a new face cycle of |f| and 
   makes |e| the entry point of it.}*/
{
  SHalfedge_around_sface_circulator hfc(e), hend(hfc);
  CGAL_For_all(hfc,hend) hfc->incident_sface() = f;
  store_sm_boundary_object(e,f);
} 

void link_as_loop(SHalfloop_handle l, SFace_handle f)
/*{\Mop creates a new trivial face cycle of |f| and 
   makes |l| the singular object of it.}*/
{ store_sm_boundary_object(l,f); l->incident_sface() = f; } 

void link_as_isolated_vertex(SVertex_handle v, SFace_handle f)
/*{\Mop creates a new trivial face cycle of |f|.
   (makes |v| an isolated vertex within |f|).}*/
{ store_sm_boundary_object(v,f); v->incident_sface() = f; } 

void unlink_as_face_cycle(SHalfedge_handle e)
/*{\Mop removes the face cycle defined by |e| from |e->incident_sface()|.
    Does not update the face links of the corresponding face cycle
    edges. \precond |e| is the entry object of the face cycle.}*/
{ undo_sm_boundary_object(e,e->incident_sface()); }
  
void unlink_as_loop(SHalfloop_handle l)
/*{\Mop removes the trivial face cycle defined by |l| from
   |l->incident_sface()|. Does not update |l|'s face link.}*/
{ undo_sm_boundary_object(l,l->incident_sface()); }

void unlink_as_isolated_vertex(SVertex_handle v)
/*{\Mop removes the trivial face cycle defined by |v| from
   |v->incident_sface()|. Does not update |v|'s face link.
   \precond |v| is a trivial face cycle of |v->incident_sface()|.}*/
{ undo_sm_boundary_object(v,v->incident_sface()); }

void clear_face_cycle_entries(SFace_handle f)
{ map()->reset_sm_object_list(f->boundary_entry_objects());
  // removes entries of list and the hashed membership
}

SHalfedge_handle new_shalfedge_pair(SVertex_handle v1,
			            SVertex_handle v2)
/*{\Mop creates a new pair of edges |(e1,e2)| representing |(v1,v2)| 
  by appending the |ei| to |A(vi)| $(i=1,2)$.}*/
{ SHalfedge_handle e1 = new_shalfedge_pair();
  SHalfedge_handle e2 = e1->twin();
  if (!is_isolated(v1))
    set_adjacency_at_source_between(cap(first_out_edge(v1)),e1,
                                    first_out_edge(v1));
  else
    close_tip_at_source(e1,v1);
  if (!is_isolated(v2))
    set_adjacency_at_source_between(cap(first_out_edge(v2)),e2,
                                    first_out_edge(v2));
  else 
    close_tip_at_source(e2,v2);
  return e1;
}


SHalfedge_handle new_shalfedge_pair(SHalfedge_handle e1, 
			            SHalfedge_handle e2,
			            int pos1 = AFTER, int pos2 = AFTER)
/*{\Mop creates a new pair of edges |(es1,es2)| representing the uedge
  |\{e1->source(),e2->source()\}| by inserting the |esi| before or after |ei| 
  into the cyclic adjacency list of |ei->source()| depending on |posi| 
  $(i=1,2)$ from |\Mname::BEFORE|, |\Mname::AFTER|.}*/
{ 
  SHalfedge_handle er = new_shalfedge_pair();
  SHalfedge_handle ero = er->twin();
  if (pos1 < 0) { // before e1
    set_adjacency_at_source_between(cap(e1),er,e1);
    if ( e1 == first_out_edge(e1->source()) )
      set_first_out_edge(e1->source(),er);
  } else { // after e1
    set_adjacency_at_source_between(e1,er,cas(e1));
  }
  if (pos2 < 0) { // before e2
    set_adjacency_at_source_between(cap(e2),ero,e2);
    if ( e2 == first_out_edge(e2->source()) )
      set_first_out_edge(e2->source(),ero);
  } else { // after e2
    set_adjacency_at_source_between(e2,ero,cas(e2));
  }
  return er;
}

SHalfedge_handle new_shalfedge_pair(SHalfedge_handle e, SVertex_handle v,
                               int pos = AFTER)
/*{\Mop creates a new pair of edges  |(e1,e2)| representing the uedge
  |\{e->source(),v\}| by inserting |e1| before or after |e| into cyclic 
  adjacency list of |e->source()| depending on |pos| from |\Mname::BEFORE|,
  |\Mname::AFTER| and appending |e2| at |A(v)|.}*/
{
  SHalfedge_handle e_new = new_shalfedge_pair();
  SHalfedge_handle e_opp = e_new->twin();
  if (pos < 0) { // before e
    set_adjacency_at_source_between(cap(e),e_new,e);
    if ( e == first_out_edge(e->source()) )
      set_first_out_edge(e->source(),e_new);
  } else  // after e
    set_adjacency_at_source_between(e,e_new,cas(e));
  
  if (!is_isolated(v)) {
    SHalfedge_handle e_first = first_out_edge(v);
    set_adjacency_at_source_between(cap(e_first),e_opp,e_first);
  } else
    close_tip_at_source(e_opp,v);
  return e_new;
}


SHalfedge_handle new_shalfedge_pair(SVertex_handle v, SHalfedge_handle e,
                                    int pos = AFTER)
/*{\Mop symmetric to the previous one.}*/
{ return new_shalfedge_pair(e,v,pos)->twin(); }

void delete_edge_pair(SHalfedge_handle e)
/*{\Mop deletes |e| and its twin and maintains the adjacency at its source 
        and its target.}*/
{ remove_from_adj_list_at_source(e);
  remove_from_adj_list_at_source(e->twin());
  delete_edge_pair_only(e);
}

SHalfedge_handle split_at(SHalfedge_handle e, Sphere_point sp)
{
  CGAL_assertion(sp != e->source()->point());
  CGAL_assertion(sp != e->twin()->source()->point());
  CGAL_assertion(Sphere_segment(e->source()->point(),
  				e->twin()->source()->point(),
  				e->circle()).has_on(sp));
  SVertex_handle v_new = new_svertex(sp);
  v_new->mark() = e->mark();
  return split_at(e, v_new);
}  

SHalfedge_handle split_at(SHalfedge_handle e, SVertex_handle v)
{
  CGAL_assertion(v->point() != e->source()->point());
  CGAL_assertion(v->point() != e->twin()->source()->point());
  CGAL_assertion(Sphere_segment(e->source()->point(),
  				e->twin()->source()->point(),
  				e->circle()).has_on(v->point()));

  SHalfedge_handle e_new = new_shalfedge_pair(v, e->twin());
  e_new->mark() = e_new->twin()->mark() = e->mark();
  e_new->circle() = e->circle();
  e_new->twin()->circle() = e->twin()->circle();
  if(e->twin()->source()->out_sedge() == e->twin())
    e->twin()->source()->out_sedge() = e_new->twin();
  e->twin()->source() = v;
  e->snext()->sprev() = e_new;
  e_new->snext() = e->snext();
  e_new->sprev() = e;
  e->snext() = e_new;
  e_new->twin()->snext() = e->twin();
  e->twin()->sprev() = e_new->twin();
  e_new->incident_sface() = e->incident_sface();
  e_new->twin()->incident_sface() = e->twin()->incident_sface();
  return e_new;
}

void delete_vertex(SVertex_handle v)
/*{\Mop deletes |v| and all outgoing edges |A(v)| as well as their twins. 
   Updates the adjacency at the targets of the edges in |A(v)|.}*/
{ 
  if (!is_isolated(v)) {
    SHalfedge_handle e = first_out_edge(v);
    while ( e != cap(e) ) 
      delete_edge_pair(cap(e));  
    delete_edge_pair(e); 
  }
  delete_vertex_only(v);
}

void delete_face(SFace_handle f)
/*{\Mop deletes the face |f| without consideration of topological 
   linkage.}*/
{ clear_face_cycle_entries(f); delete_face_only(f); }

bool has_outdeg_two(SVertex_handle v) const
/*{\Mop return true when |v| has outdegree two.}*/
// does this work for looping edges?
{ if (is_isolated(v)) return false;
  SHalfedge_handle e1 = first_out_edge(v);
  SHalfedge_handle e2 = last_out_edge(v);
  return (e1!=e2 && e2==cas(e1));
}

void link_as_prev_next_pair(SHalfedge_handle e1, SHalfedge_handle e2)
/*{\Xop makes |e1| and |e2| adjacent in the face cycle 
   $\ldots-|e1-e2|-\ldots$.
   Afterwards |e1 = e2->sprev()| and |e2 = e1->snext()|.}*/
{ e1->snext() = e2; e2->sprev() = e1; }

void merge_edge_pairs_at_target(SHalfedge_handle e)
/*{\Mop merges the edge pairs at |v = e->target()|. |e| and |twin(e)| 
  are preserved, |e->snext()|, |twin(e->snext())| and |v| are deleted
  in the merger. \precond |v| has outdegree two. The adjacency at 
  |e->source()| and |target(e->snext())| is kept consistent.
  If |e->snext()| was entry point of |e->incident_sface()| then |e| takes this role.
  The same holds for |twin(e->snext())| and |face(twin(e))|.}*/
{
  CGAL_NEF_TRACEN("merge_edge_pairs_at_target "<<PH(e));
  SHalfedge_handle en = e->snext(), eno = en->twin(), enn, enno,
               eo = e->twin() ;
  if ( is_closed_at_target(en) ) { enn = eo; enno=e; }
  else { enn = en->snext(), enno = eno->sprev(); }
  SVertex_handle v = e->target(), vn = en->target();
  CGAL_assertion(has_outdeg_two(v));
  SFace_handle f1 = en->incident_sface(), f2 = eno->incident_sface();
  // transfer the opposite face cycles e-en-enn to e-enn
  if ( enn != eno ) {
    link_as_prev_next_pair(e,enn);
    link_as_prev_next_pair(enno,eo);
  } else {
    link_as_prev_next_pair(e,eo);
  }
  // set vertex of e and deal with vertex-halfedge incidence
  eo->source() = vn;

  if ( first_out_edge(vn) == eno ) set_first_out_edge(vn,eo);
  if ( is_sm_boundary_object(en) )
  { undo_sm_boundary_object(en,f1); store_sm_boundary_object(e,f1); }
  if ( is_sm_boundary_object(eno) )
  { undo_sm_boundary_object(eno,f2); store_sm_boundary_object(eo,f2); }
  delete_vertex_only(v);
  delete_edge_pair_only(en);
  CGAL_NEF_TRACEN("END "<<PH(e->sprev())<<PH(e)<<PH(e->snext()));
}

void convert_edge_to_loop(SHalfedge_handle e)
/*{\Mop converts the edge at |v = e->target()| to the unique
  loop |l| of |\Mvar|. |e|, |e->twin()| and |v| are deleted
  in the conversion. \precond there was no loop in |\Mvar|.
  As |e| was entry point of |e->incident_sface()| then |l| takes this role.}*/
{ CGAL_NEF_TRACEN("convert_edge_to_loop "<<PH(e));
  CGAL_assertion( e->source()==e->target() );
  CGAL_assertion( !has_shalfloop() );
  SHalfloop_handle l = new_shalfloop_pair();
  SVertex_handle v = e->target();
  SFace_handle f1 = e->incident_sface(), f2 = e->twin()->incident_sface();
  if( is_sm_boundary_object(e)) {
    CGAL_assertion( is_sm_boundary_object(e->twin()));
    undo_sm_boundary_object(e,f1); undo_sm_boundary_object(e->twin(),f2);
  }
  link_as_loop(l,f1), link_as_loop(l->twin(),f2);
  l->circle() = e->circle(); l->twin()->circle() = e->twin()->circle();
  l->mark() = l->twin()->mark() = e->mark();
  delete_vertex_only(v);
  delete_edge_pair_only(e);
}

void flip_diagonal(SHalfedge_handle e)
{ SHalfedge_handle r = e->twin();
  SHalfedge_handle en = e->snext(), enn= en->snext();
  SHalfedge_handle rn = r->snext(), rnn= rn->snext();
  CGAL_NEF_TRACEN(PH(e)<<PH(en)<<PH(enn));
  CGAL_NEF_TRACEN(PH(r)<<PH(rn)<<PH(rnn));
  CGAL_assertion( enn->snext()==e && rnn->snext()==r );
  remove_from_adj_list_at_source(e);
  remove_from_adj_list_at_source(r);
  set_adjacency_at_source_between(enn,e,en->twin());
  set_adjacency_at_source_between(rnn,r,rn->twin());
}
     

/*{\Mtext \headerline{Incomplete topological update primitives}}*/

SHalfedge_handle new_shalfedge_pair_at_source
  (SHalfedge_handle e, int pos = AFTER) const
/*{\Xop creates a new pair of edges  |(e1,e2)| representing |(e->source(),())| 
  by inserting |e1| before or after |e| into cyclic adjacency list of
  |e->source()| depending on |pos| from |\Mname::BEFORE, \Mname::AFTER|.}*/
{ SHalfedge_handle e_new = new_shalfedge_pair();
  if (pos < 0) set_adjacency_at_source_between(cap(e),e_new,e);
  else         set_adjacency_at_source_between(e,e_new,cas(e));
  return e_new;
}

SHalfedge_handle new_shalfedge_pair_at_source
  (SVertex_handle v, int pos = AFTER)
/*{\Mop creates a new pair of edges  |(e1,e2)| representing |(v,())| 
  by inserting |e1| at the beginning (BEFORE) or end (AFTER)
  of adjacency list of |v|.}*/
{ SHalfedge_handle e1 = new_shalfedge_pair();
  if (!is_isolated(v)) {
    SHalfedge_handle ef = first_out_edge(v);
    if (pos < 0) { // before e1
      set_adjacency_at_source_between(cap(ef),e1,ef);
      set_first_out_edge(v,e1);
    } else {
      set_adjacency_at_source_between(cap(ef),e1,ef);
    }
  } else
    close_tip_at_source(e1,v);
  return e1;
}

void delete_edge_pair_at_source(SHalfedge_handle e)
/*{\Mop deletes |e| and its twin and maintains the adjacency at its 
  source.}*/
{ remove_from_adj_list_at_source(e);
  delete_edge_pair_only(e);
}

void link_as_target_and_append(SVertex_handle v, SHalfedge_handle e, int pos = AFTER)
/*{\Mop makes |v| the target of |e| appends |e->twin()| to the adjacency list
   of |v|.}*/
{ if (!is_isolated(v)) {
    SHalfedge_handle ef = first_out_edge(v);
    set_adjacency_at_source_between(cap(ef), e->twin(), ef);
    if(pos<0)
      set_first_out_edge(v,e->twin());
  } else
    close_tip_at_target(e,v);
}

void link_as_source_of(SHalfedge_handle e, SVertex_handle v) const
/*{\Mop makes |e->source() = v| and sets |e| as the first
        out edge if |v| was isolated before.}*/
{ e->source() = v;
  if (v->out_sedge() == SHalfedge_handle()) v->out_sedge() = e; }

void link_as_target_of(SHalfedge_handle e, SVertex_handle v) const
/*{\Mop makes |e->target() = v| and sets |e| as the first
  in edge if |v| was isolated before.}*/
{ link_as_source_of(e->twin(),v); }

void set_adjacency_at_source_between(SHalfedge_handle e, SHalfedge_handle en)
/*{\Mop makes |e| and |en| neigbors in the cyclic ordered adjacency list 
    around |v=e->source()|. \precond |e->source()==en->source()|.}*/
{ CGAL_assertion(e->source()==en->source());
  link_as_prev_next_pair(en->twin(),e);
}

void set_adjacency_at_source_between(SHalfedge_handle e1, 
                                     SHalfedge_handle e_between, 
                                     SHalfedge_handle e2)
/*{\Mop inserts |e_between| into the adjacency list around |e1->source()| 
  between |e1| and |e2| and makes |e1->source()| the source of |e_between|. 
  \precond |e1->source()==e2->source()|.}*/
{ 
  e_between->source() = e1->source();
  set_adjacency_at_source_between(e1,e_between);
  set_adjacency_at_source_between(e_between,e2);
}

void close_tip_at_source(SHalfedge_handle e, SVertex_handle v) 
/*{\Mop sets |v| as source of |e| and closes the tip by setting the 
  corresponding pointers such that |prev(e) == e->twin()| and
  |next(e->twin()) == e|.}*/
{ link_as_source_of(e,v); 
  link_as_prev_next_pair(e->twin(),e); }

void close_tip_at_target(SHalfedge_handle e, SVertex_handle v) 
/*{\Mop sets |v| as target of |e| and closes the tip by setting the 
  corresponding pointers such that |prev(e->twin()) == e| and 
  |e->snext() == e->twin()|.}*/
{ link_as_target_of(e,v);
  link_as_prev_next_pair(e,e->twin()); }


void remove_from_adj_list_at_source(SHalfedge_handle e)
/*{\Mop removes a halfedge pair |(e,e->twin()| from the adjacency list
of |e->source()|. Afterwards |next(prev(e))==next(e->twin())| and
|first_out_edge(e->source())| is valid if |degree(v->source())>1| before
the operation.}*/
{
  SVertex_handle v = e->source();
  if ( is_closed_at_source(e) ) { // last outgoing
    v->out_sedge() = SHalfedge_handle();
  } else {
    if (e == first_out_edge(v)) v->out_sedge() = cap(e);
    set_adjacency_at_source_between(cap(e),cas(e));
  }
}

void set_face(SHalfedge_handle e, SFace_handle f) const
{ e->incident_sface() = f; }
void set_face(SHalfloop_handle l, SFace_handle f) const
{ l->incident_sface() = f; }
void set_face(SVertex_handle v, SFace_handle f) const
{ v->incident_sface() = f; }
void set_first_out_edge(SVertex_handle v, SHalfedge_handle e) const
{ v->out_sedge() = e; }
void set_prev(SHalfedge_handle e, SHalfedge_handle ep) const
{ e->sprev() = ep; }
void set_next(SHalfedge_handle e, SHalfedge_handle en) const
{ e->snext() = en; }
void set_source(SHalfedge_handle e, SVertex_handle v) const
{ e->source() = v; }

/*{\Mtext \headerline{Associated Information}\restoreopdims}*/

void set_marks_in_face_cycle(SHalfedge_handle e, Mark m) const
{ SHalfedge_around_sface_circulator hfc(e), hend(hfc);
  CGAL_For_all(hfc,hend) hfc->mark() = hfc->target()->mark() = m;
}

GenPtr& info(SVertex_handle v) const
{ return v->info(); }
GenPtr& info(SHalfedge_handle e) const
{ return e->info(); }
GenPtr& info(SFace_handle f) const
{ return f->info(); }

const GenPtr& info(SVertex_const_handle v) const
{ return v->info(); }
const GenPtr& info(SHalfedge_const_handle e) const
{ return e->info(); }
const GenPtr& info(SFace_const_handle f) const
{ return f->info(); }

 
/*{\Mtext \headerline{Iteration}}*/
/*{\Mtext The list of all objects can be accessed via iterator ranges.
For comfortable iteration we also provide iterations macros. 
The iterator range access operations are of the following kind:\\
|SVertex_iterator svertices_begin()/svertices_end()|\\
|SHalfedge_iterator shalfedges_begin()/shalfedges_end()|\\
|SHalfedge_iterator sedges_begin()/sedges_end()|\\
|SFace_iterator sfaces_begin()/sfaces_end()|

The macros are then |CGAL_forall_svertices_of(v,V)|,
|CGAL_forall_shalfedges_of(e,V)|, |CGAL_forall_sedges_of(e,V)|, 
|CGAL_forall_sfaces_of(f,V)|, |CGAL_forall_sface_cycles_of(fc,F)|.}*/

void transform( const Aff_transformation_3& linear) {
  //  CGAL_NEF_TRACEN("transform sphere map of vertex" << center_vertex()->point());
    // The affine transformation is linear, i.e., no translation part.
    CGAL_precondition( linear.hm(0,3) == 0 && 
                       linear.hm(1,3) == 0 && 
                       linear.hm(2,3) == 0);

    //    CGAL_NEF_TRACEN(linear);
    for (SVertex_iterator i = svertices_begin(); i != svertices_end(); ++i)
      i->point() = normalized(Sphere_point( i->point().transform( linear)));
    for (SHalfedge_iterator i = shalfedges_begin(); i !=shalfedges_end(); ++i)
      i->circle() = Sphere_circle( i->circle().transform( linear));
    if ( has_shalfloop()) {
      shalfloop()->circle() = Sphere_circle(shalfloop()->circle()
					  .transform( linear));
      shalfloop()->twin()->circle()
	= Sphere_circle(shalfloop()->twin()->circle().transform( linear));
    }
}

void extract_complement() {

  SVertex_handle sv;
  CGAL_forall_svertices(sv,*this) sv->mark() = !sv->mark();
  SHalfedge_handle she;
  CGAL_forall_shalfedges(she,*this) she->mark() = !she->mark();
  SHalfloop_handle shl;
  if(has_shalfloop()) {
    shl = shalfloop(); 
    shl->mark() = shl->twin()->mark() = !shl->mark();
  }
  SFace_handle sf;
  CGAL_forall_sfaces(sf,*this) 
    sf->mark() = !sf->mark();
}

void extract_interior() {

  SVertex_handle sv;
  CGAL_forall_svertices(sv,*this) sv->mark() = false;
  SHalfedge_handle she;
  CGAL_forall_shalfedges(she,*this) she->mark() = false;
  SHalfloop_handle shl;
  if(has_shalfloop()) { 
    shl = shalfloop(); 
    shl->mark() = shl->twin()->mark() = false;
  }
}

void extract_boundary() {

  SVertex_handle sv;
  CGAL_forall_svertices(sv,*this) sv->mark() = true;
  SHalfedge_handle she;
  CGAL_forall_shalfedges(she,*this) she->mark() = true;
  SHalfloop_handle shl;
  if(has_shalfloop()) { 
    shl = shalfloop(); 
    shl->mark() = shl->twin()->mark() = true;
  }
  SFace_handle sf;
  CGAL_forall_sfaces(sf,*this) sf->mark() = false;
}

bool is_valid( Unique_hash_map<SVertex_handle,bool>& sv_hash,
	       Unique_hash_map<SHalfedge_handle,bool>& se_hash,
	       Unique_hash_map<SFace_handle,bool>& sf_hash,
	       bool verb = false, int level = 0) {
    
  Verbose_ostream verr(verb);
  verr << "begin CGAL::SNC_SM_decorator<...>::is_valid( verb=true, "
    "level = " << level << "):" << std::endl;
    
  bool valid = true;

  std::size_t count = 0;
  std::size_t max = 2 * number_of_svertices() 
    + 2 * number_of_shalfedges()
    + number_of_sfaces()
    + 2;
  
  SVertex_handle sv;
  int isolated_vertices_found = 0;
  CGAL_forall_svertices(sv,*this) {
    valid = valid && (!sv_hash[sv]);
    sv_hash[sv] = true;
    if(is_isolated(sv)) 
      isolated_vertices_found++;
    valid = valid && (++count <= max);
  }

  SHalfedge_iterator she;
  CGAL_forall_shalfedges(she,*this) {
    valid = valid && she->is_valid(verb, level);

    valid = valid && (she->twin() != she);
    valid = valid && (she->twin()->twin() == she);
    valid = valid && (she->snext()->sprev() == she);
    valid = valid && ((she->sprev() != she && she->snext() != she) || 
		      (she->sprev() == she && she->snext() == she));
    valid = valid && (she->incident_sface() == she->snext()->incident_sface());
    valid = valid && (she->incident_sface() == she->sprev()->incident_sface());

    valid = valid && (!se_hash[she]);

    //    Plane_3 pl(point(she->source()),point(she->target()),Point_3(0,0,0));
    //    Sphere_point vct(pl.orthogonal_vector());
    //    valid = valid && (normalized(Sphere_point(she->circle().orthogonal_vector())) == normalized(vct) || 
    //	      normalized(Sphere_point(she->circle().opposite().orthogonal_vector())) == normalized(vct));

    se_hash[she] = true;
    valid = valid && (++count <= max);
  }

  if(has_shalfloop()) {
    SHalfloop_handle shl = shalfloop();
    valid = valid && shl->is_valid();
    valid = valid && shl->twin()->is_valid();
    valid = valid && (shl->twin() != shl);
    valid = valid && (shl->twin()->twin() == shl); 
  }

  SFace_iterator sf;
  SFace_cycle_iterator sfc;  
  int loop_entries_found = 0;
  int edge_entries_found = 0;
  int vertex_entries_found = 0;
  CGAL_forall_sfaces(sf,*this) {
    valid = valid && sf->is_valid(verb, level);
    valid = valid && (!sf_hash[sf]);
    sf_hash[sf] = true;
    
    CGAL_forall_sface_cycles_of(sfc,sf) {
      if (sfc.is_shalfloop()) 
	loop_entries_found++;
      else if(sfc.is_shalfedge())
	edge_entries_found++;
      else if(sfc.is_svertex())
	vertex_entries_found++;
      valid = valid && (++count <= max);
    }
    
    valid = valid && (++count <= max);
  }

  if(has_shalfloop())
    valid = valid && (loop_entries_found == 2);
  else
    valid = valid && (loop_entries_found == 0);

  if(number_of_shalfedges() > 0)
    valid = valid && (edge_entries_found > 0);
  else
    valid = valid && (edge_entries_found == 0);

  valid = valid && (vertex_entries_found == isolated_vertices_found);

  verr << "end of CGAL::SNC_SM_decorator<...>::is_valid(): structure is "
       << ( valid ? "valid." : "NOT VALID.") << std::endl;

  return valid;
}

  template <typename Selection>
  void change_marks(const Mark& m, const Selection& SP) {
   
    psm_->mark() = SP(m, psm_->mark());

    SVertex_iterator v;
    CGAL_forall_svertices(v,*this)
      v->mark() = SP(m, v->mark());

    SHalfedge_iterator e;
    CGAL_forall_shalfedges(e,*this)
      e->mark() = SP(m, e->mark());

    SFace_iterator f;
    CGAL_forall_sfaces(f,*this) 
      f->mark() = SP(m, f->mark());
  }

  template <typename Selection>
  void change_marks(const Selection& SP, const Mark& m) {
   
    psm_->mark() = SP(psm_->mark(), m);

    SVertex_iterator v;
    CGAL_forall_svertices(v,*this)
      v->mark() = SP(v->mark(),m);

    SHalfedge_iterator e;
    CGAL_forall_shalfedges(e,*this)
      e->mark() = SP(e->mark(),m);

    SFace_iterator f;
    CGAL_forall_sfaces(f,*this) 
      f->mark() = SP(f->mark(),m);
  }

void check_integrity_and_topological_planarity(bool faces=true) const {
  SM_const_decorator C(psm_);
  C.check_integrity_and_topological_planarity(faces);
}

}; // SM_decorator


} //namespace CGAL
#endif // CGAL_SM_DECORATOR_H
