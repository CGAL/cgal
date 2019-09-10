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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_SNC_SM_OVERLAYER_H
#define CGAL_SNC_SM_OVERLAYER_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>
#include <CGAL/Union_find.h>
#include <CGAL/Nef_2/Segment_overlay_traits.h>
#ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <CGAL/Nef_2/geninfo.h>
#else
#include <boost/any.hpp>
#endif
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_const_decorator.h>
#include <CGAL/Nef_S2/SM_overlayer.h>
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 131
#include <CGAL/Nef_2/debug.h>

#ifndef CGAL_USE_LEDA
#define LEDA_MEMORY(t) 
#endif

namespace CGAL {

/*{\Manpage {SNC_SM_overlayer}{Refs_}{Overlay in the sphere}{O}}*/

template <typename Items, typename SM_decorator_>
class SNC_SM_overlayer : public SM_overlayer<SM_decorator_> {
public:

  typedef SM_decorator_                          SM_decorator;
  typedef typename SM_decorator::Map             Map;
  typedef SM_overlayer<SM_decorator_>            Base;
  typedef SNC_SM_overlayer<Items, SM_decorator_> Self;

  typedef typename Base::SVertex_handle SVertex_handle;
  typedef typename Base::SHalfedge_handle SHalfedge_handle;
  typedef typename Base::SHalfloop_handle SHalfloop_handle;
  typedef typename Base::SFace_handle SFace_handle;
  typedef typename Base::SVertex_iterator SVertex_iterator;
  typedef typename Base::SHalfedge_iterator SHalfedge_iterator;
  typedef typename Base::SFace_iterator SFace_iterator;
  typedef typename Base::SHalfedge_around_sface_circulator 
                         SHalfedge_around_sface_circulator;

  typedef typename Base::Sphere_kernel           Sphere_kernel;

  typedef typename Map::Infi_box Infi_box;

  using Base::clear_face_cycle_entries;
  using Base::link_as_loop;
  using Base::is_closed_at_source;
  using Base::is_isolated;
  using Base::delete_edge_pair;
  using Base::delete_vertex_only;
  using Base::delete_face_only;
  using Base::store_sm_boundary_object;
  using Base::first_out_edge;
  using Base::set_face;
  using Base::has_outdeg_two;
  using Base::convert_edge_to_loop;
  using Base::merge_edge_pairs_at_target;

 public:
  SNC_SM_overlayer(Map* M, 
		   const Sphere_kernel& G = Sphere_kernel()) 
    : Base(M,G) {}

  template <typename Association>
  void simplify(Association&) {
    CGAL_NEF_TRACEN("simplifying"); 
    
    typedef typename CGAL::Union_find<SFace_handle>::handle Union_find_handle;
    CGAL::Unique_hash_map< SFace_handle, Union_find_handle> Pitem(NULL);
    CGAL::Unique_hash_map< SVertex_handle, Union_find_handle> Vitem(NULL);
    CGAL::Union_find< SFace_handle> UF;
    
    SFace_iterator f;
    CGAL_forall_sfaces(f,*this) {
      Pitem[f] = UF.make_set(f);
      clear_face_cycle_entries(f);
    }
    
    if ( this->has_shalfloop() ) {
      SHalfloop_handle l = this->shalfloop();
      SFace_handle f = *(UF.find(Pitem[l->incident_sface()]));
      link_as_loop(l,f);
      f = *(UF.find(Pitem[l->twin()->incident_sface()]));
      link_as_loop(l->twin(),f);
    }
    
    SHalfedge_iterator e, en;
    for(e = this->shalfedges_begin(); e != this->shalfedges_end(); e = en) { 
      en = e; ++en; if ( en==e->twin() ) ++en;
      CGAL_NEF_TRACEN("can simplify ? " << PH(e));
      if(!Infi_box::is_sedge_on_infibox(e)) {
	CGAL_NEF_TRACEN(e->mark() << " " << 
			e->incident_sface()->mark() << " " << 
			e->twin()->incident_sface()->mark());
	if (( e->mark() == e->incident_sface()->mark() && 
	      e->incident_sface()->mark() == e->twin()->incident_sface()->mark())){
	  CGAL_NEF_TRACEN("deleting "<<PH(e));
	  if ( !UF.same_set(Pitem[e->incident_sface()],
			    Pitem[e->twin()->incident_sface()]) ) {
	    
	    UF.unify_sets( Pitem[e->incident_sface()],
			   Pitem[e->twin()->incident_sface()] );
	    CGAL_NEF_TRACEN("unioning disjoint faces");
	  }
	  
	  CGAL_NEF_TRACEN("is_closed_at_source " << is_closed_at_source(e) << 
			  " " << is_closed_at_source(e->twin()));
	  
	  if ( is_closed_at_source(e) )
	    Vitem[e->source()] = Pitem[e->incident_sface()];
	  
	  if ( is_closed_at_source(e->twin()))
	    Vitem[e->target()] = Pitem[e->incident_sface()];
	  
	  delete_edge_pair(e);
	}
      }
    }
    
    CGAL::Unique_hash_map<SHalfedge_handle,bool> linked(false);
    for (e = this->shalfedges_begin(); e != this->shalfedges_end(); ++e) {
      if ( linked[e] ) continue;
      SHalfedge_around_sface_circulator hfc(e),hend(hfc);
      SFace_handle f = *(UF.find( Pitem[e->incident_sface()]));
      CGAL_For_all(hfc,hend) {  set_face(hfc,f); linked[hfc]=true; }
      store_sm_boundary_object(e,f);
    }
    
    SVertex_iterator v,vn;
    for(v = this->svertices_begin(); v != this->svertices_end(); v=vn) {
      vn=v; ++vn;
      if ( is_isolated(v) ) {
	if(Vitem[v] != NULL) {
	  set_face(v,*(UF.find(Vitem[v])));
	  CGAL_NEF_TRACEN("incident face of " << PH(v) << " set to " << &*(v->incident_sface()));
	}
	else {
	set_face(v, *(UF.find(Pitem[v->incident_sface()])));
	CGAL_NEF_TRACEN("isolated svertex " << PH(v) << 
			" already has incident face " << &*(v->incident_sface()));
	}
	if ( v->mark() == v->incident_sface()->mark() ) {
	  CGAL_NEF_TRACEN("removing isolated vertex"<<PH(v));
	  delete_vertex_only(v);  
	} 
	else 
	  store_sm_boundary_object(v,v->incident_sface()); // isolated, but should stay
      } else { // v not isolated
	SHalfedge_handle e2 = first_out_edge(v), e1 = e2->sprev();
	if ( has_outdeg_two(v) &&
	     v->mark() == e1->mark() && e1->mark() == e2->mark() &&
	     e1->circle() == e2->circle() ) {
	  CGAL_NEF_TRACEN("collinear at "<<PH(v)<<PH(e1)<<PH(e2));
	  if ( e1 == e2 ){ 
	    CGAL_NEF_TRACEN("edge_to_loop"); 
	    convert_edge_to_loop(e1);
	  } else {
	    CGAL_NEF_TRACEN("merge_edge_pairs"); 
	    merge_edge_pairs_at_target(e1); 
	  } 	
	}
      }
    }
    
    SFace_iterator fn;
    for (f = fn = this->sfaces_begin(); f != this->sfaces_end(); f=fn) { 
      ++fn;
      Union_find_handle pit = Pitem[f];
      if ( UF.find(pit) != pit ) {
	CGAL_NEF_TRACEN("delete face " << &*f);
	delete_face_only(f);
      }
    }
    
    CGAL_NEF_TRACEN(" ");
    CGAL_NEF_TRACEN("resulting vertex ");
    
    CGAL_forall_svertices(vn, *this)
      CGAL_NEF_TRACEN("|" << vn->point() << "|" << vn->mark());
    CGAL_NEF_TRACEN(" ");
    
    CGAL_forall_shalfedges(en,*this)
      CGAL_NEF_TRACEN("|" << en->circle() <<
		      "|" << en->mark() << 
		      " " << en->incident_sface()->mark());
    CGAL_NEF_TRACEN("---------------------");
  }

};

template <typename SM_decorator_>
class SNC_SM_overlayer<SNC_indexed_items, SM_decorator_> 
  : public SM_overlayer<SM_decorator_> {

  typedef SM_decorator_                          SM_decorator;
  typedef typename SM_decorator::Map             Map;
  typedef SM_overlayer<SM_decorator_>            Base;
  typedef SNC_SM_overlayer<SNC_indexed_items, 
    SM_decorator_> Self;

  typedef typename Base::SVertex_handle SVertex_handle;
  typedef typename Base::SHalfedge_handle SHalfedge_handle;
  typedef typename Base::SHalfloop_handle SHalfloop_handle;
  typedef typename Base::SFace_handle SFace_handle;
  typedef typename Base::SVertex_iterator SVertex_iterator;
  typedef typename Base::SHalfedge_iterator SHalfedge_iterator;
  typedef typename Base::SFace_iterator SFace_iterator;
  typedef typename Base::SHalfedge_around_sface_circulator 
                         SHalfedge_around_sface_circulator;

  typedef typename Base::Sphere_kernel           Sphere_kernel;

  typedef typename Map::Infi_box Infi_box;

  using SM_decorator::clear_face_cycle_entries;
  using SM_decorator::link_as_loop;
  using SM_decorator::link_as_prev_next_pair;
  using SM_decorator::is_closed_at_source;
  using SM_decorator::is_closed_at_target;
  using SM_decorator::delete_edge_pair;
  using SM_decorator::set_face;
  using SM_decorator::is_isolated;
  using SM_decorator::delete_vertex_only;
  using SM_decorator::delete_face_only;
  using SM_decorator::first_out_edge;
  using SM_decorator::set_first_out_edge;
  using SM_decorator::has_outdeg_two;
  using SM_decorator::store_sm_boundary_object;
  using SM_decorator::is_sm_boundary_object;
  using SM_decorator::undo_sm_boundary_object;
  using SM_decorator::delete_edge_pair_only;

 public:
  SNC_SM_overlayer(Map* M, 
		   const Sphere_kernel& G = Sphere_kernel()) 
    : Base(M,G) {}

  void convert_edge_to_loop(SHalfedge_handle e) {
    /*{\Mop converts the edge at |v = e->target()| to the unique
      loop |l| of |\Mvar|. |e|, |e->twin()| and |v| are deleted
      in the conversion. \precond there was no loop in |\Mvar|.
      As |e| was entry point of |e->incident_sface()| then |l| takes this role.}*/
    
    CGAL_NEF_TRACEN("convert_edge_to_loop "<<PH(e));
    CGAL_assertion( e->source()==e->target() );
    CGAL_assertion( !this->has_shalfloop() );
    SHalfloop_handle l = this->new_shalfloop_pair();
    SVertex_handle v = e->target();
    SFace_handle f1 = e->incident_sface(), f2 = e->twin()->incident_sface();
    if( is_sm_boundary_object(e)) {
      CGAL_assertion( is_sm_boundary_object(e->twin()));
      undo_sm_boundary_object(e,f1); undo_sm_boundary_object(e->twin(),f2);
    }
    link_as_loop(l,f1), link_as_loop(l->twin(),f2);
    l->circle() = e->circle(); l->twin()->circle() = e->twin()->circle();
    l->mark() = l->twin()->mark() = e->mark();
    l->set_index(e->get_index());
    l->twin()->set_index(e->twin()->get_index());
    delete_vertex_only(v);
    delete_edge_pair_only(e);
  }
  
  template <typename Association>
  void merge_edge_pairs_at_target(SHalfedge_handle e, Association& A) {
    /*{\Mop merges the edge pairs at |v = e->target()|. |e| and |twin(e)| 
      are preserved, |e->snext()|, |twin(e->snext())| and |v| are deleted
      in the merger. \precond |v| has outdegree two. The adjacency at 
      |e->source()| and |target(e->snext())| is kept consistent.
      If |e->snext()| was entry point of |e->incident_sface()| then |e| takes this role.
      The same holds for |twin(e->snext())| and |face(twin(e))|.}*/

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
    
    CGAL_NEF_TRACEN("rehash " << en->get_index() << " " << e->get_index());
    CGAL_NEF_TRACEN("       " << A.get_hash(en->get_index()) << " " << A.get_hash(e->get_index()));
    CGAL_NEF_TRACEN("rehash " << eno->get_index() << " " << eo->get_index());
    CGAL_NEF_TRACEN("       " << A.get_hash(eno->get_index()) << " " << A.get_hash(eo->get_index()));
    
    int index1 = A.get_hash(e->get_index());
    int index2 = A.get_hash(en->get_index());
    if(index2 < index1) {
      A.set_hash(e->get_index(), index2);
      e->set_index(index2);
    } else
      A.set_hash(en->get_index(), index1);

    index1 = A.get_hash(eo->get_index());
    index2 = A.get_hash(eno->get_index());    
    if(index2 < index1) {
      A.set_hash(eo->get_index(), index2);
      eo->set_index(index2);
    } else
      A.set_hash(eno->get_index(), index1);

    CGAL_NEF_TRACEN("hash sedge " << e->get_index() 
		    << "->" << A.get_hash(e->get_index()));
    CGAL_NEF_TRACEN("hash sedge " << en->get_index() 
		    << "->" << A.get_hash(en->get_index()));
    CGAL_NEF_TRACEN("hash sedge " << eo->get_index() 
		    << "->" << A.get_hash(eo->get_index()));
    CGAL_NEF_TRACEN("hash sedge " << eno->get_index() 
		    << "->" << A.get_hash(eno->get_index()));


    if ( first_out_edge(vn) == eno ) set_first_out_edge(vn,eo);
    if ( is_sm_boundary_object(en) )
      { undo_sm_boundary_object(en,f1); store_sm_boundary_object(e,f1); }
    if ( is_sm_boundary_object(eno) )
      { undo_sm_boundary_object(eno,f2); store_sm_boundary_object(eo,f2); }
    delete_vertex_only(v);
    delete_edge_pair_only(en);
    CGAL_NEF_TRACEN("END "<<PH(e->sprev())<<PH(e)<<PH(e->snext()));
  }

  template <typename Association>
  void simplify(Association& A) {
    CGAL_NEF_TRACEN("simplifying"); 
    
    typedef typename CGAL::Union_find<SFace_handle>::handle Union_find_handle;
    CGAL::Unique_hash_map< SFace_handle, Union_find_handle> Pitem(NULL);
    CGAL::Unique_hash_map< SVertex_handle, Union_find_handle> Vitem(NULL);
    CGAL::Union_find< SFace_handle> UF;
    
    SFace_iterator f;
    CGAL_forall_sfaces(f,*this) {
      Pitem[f] = UF.make_set(f);
      clear_face_cycle_entries(f);
    }
    
    if ( this->has_shalfloop() ) {
      SHalfloop_handle l = this->shalfloop();
      SFace_handle f = *(UF.find(Pitem[l->incident_sface()]));
      link_as_loop(l,f);
      f = *(UF.find(Pitem[l->twin()->incident_sface()]));
      link_as_loop(l->twin(),f);
    }
    
    SHalfedge_iterator e, en;
    for(e = this->shalfedges_begin(); e != this->shalfedges_end(); e = en) { 
      en = e; ++en; if ( en==e->twin() ) ++en;
      CGAL_NEF_TRACEN("can simplify ? " << PH(e));
      if(!Infi_box::is_sedge_on_infibox(e)) {
	CGAL_NEF_TRACEN(e->mark() << " " << 
			e->incident_sface()->mark() << " " << 
			e->twin()->incident_sface()->mark());
	if (( e->mark() == e->incident_sface()->mark() && 
	      e->incident_sface()->mark() == e->twin()->incident_sface()->mark())){
	  CGAL_NEF_TRACEN("deleting "<<PH(e));
	  if ( !UF.same_set(Pitem[e->incident_sface()],
			    Pitem[e->twin()->incident_sface()]) ) {
	    
	    UF.unify_sets( Pitem[e->incident_sface()],
			   Pitem[e->twin()->incident_sface()] );
	    CGAL_NEF_TRACEN("unioning disjoint faces");
	  }
	  
	  CGAL_NEF_TRACEN("is_closed_at_source " << is_closed_at_source(e) << 
			  " " << is_closed_at_source(e->twin()));
	  
	  if ( is_closed_at_source(e) )
	    Vitem[e->source()] = Pitem[e->incident_sface()];
	  
	  if ( is_closed_at_source(e->twin()))
	    Vitem[e->target()] = Pitem[e->incident_sface()];
	  
	  delete_edge_pair(e);
	}
      }
    }
    
    CGAL::Unique_hash_map<SHalfedge_handle,bool> linked(false);
    for (e = this->shalfedges_begin(); e != this->shalfedges_end(); ++e) {
      if ( linked[e] ) continue;
      SHalfedge_around_sface_circulator hfc(e),hend(hfc);
      SFace_handle f = *(UF.find( Pitem[e->incident_sface()]));
      CGAL_For_all(hfc,hend) {  set_face(hfc,f); linked[hfc]=true; }
      store_sm_boundary_object(e,f);
    }
    
    SVertex_iterator v,vn;
    for(v = this->svertices_begin(); v != this->svertices_end(); v=vn) {
      vn=v; ++vn;
      if ( is_isolated(v) ) {
	if(Vitem[v] != NULL) {
	  set_face(v,*(UF.find(Vitem[v])));
	  CGAL_NEF_TRACEN("incident face of " << PH(v) << " set to " << &*(v->incident_sface()));
	}
	else {
	set_face(v, *(UF.find(Pitem[v->incident_sface()])));
	CGAL_NEF_TRACEN("isolated svertex " << PH(v) << 
			" already has incident face " << &*(v->incident_sface()));
	}
	if ( v->mark() == v->incident_sface()->mark() ) {
	  CGAL_NEF_TRACEN("removing isolated vertex"<<PH(v));
	  delete_vertex_only(v);  
	} 
	else 
	  store_sm_boundary_object(v,v->incident_sface()); // isolated, but should stay
      } else { // v not isolated
	SHalfedge_handle e2 = first_out_edge(v), e1 = e2->sprev();
	if ( has_outdeg_two(v) &&
	     v->mark() == e1->mark() && e1->mark() == e2->mark() &&
	     e1->circle() == e2->circle() ) {
	  CGAL_NEF_TRACEN("collinear at "<<PH(v)<<PH(e1)<<PH(e2));
	  if ( e1 == e2 ){ 
	    CGAL_NEF_TRACEN("edge_to_loop"); 
	    convert_edge_to_loop(e1);
	  } else {
	    CGAL_NEF_TRACEN("merge_edge_pairs"); 
	    merge_edge_pairs_at_target(e1, A); 
	  } 	
	}
      }
    }
    
    SFace_iterator fn;
    for (f = fn = this->sfaces_begin(); f != this->sfaces_end(); f=fn) { 
      ++fn;
      Union_find_handle pit = Pitem[f];
      if ( UF.find(pit) != pit ) {
	CGAL_NEF_TRACEN("delete face " << &*f);
	delete_face_only(f);
      }
    }
    
    CGAL_NEF_TRACEN(" ");
    CGAL_NEF_TRACEN("resulting vertex ");
    
    CGAL_forall_svertices(vn, *this)
      CGAL_NEF_TRACEN("|" << vn->point() << "|" << vn->mark());
    CGAL_NEF_TRACEN(" ");
    
    CGAL_assertion_code(CGAL_forall_shalfedges(en,*this)
			CGAL_NEF_TRACEN("|" << en->circle() <<
					"|" << en->mark() << 
					" " << en->incident_sface()->mark()));

    CGAL_NEF_TRACEN("check indexes");
    CGAL_assertion_code(CGAL_forall_sedges(en, *this)
			CGAL_NEF_TRACEN(en->source()->point() << "->"
					<< en->twin()->source()->point() << " : "
					<< en->get_index() << " " 
					<< en->twin()->get_index()));
    CGAL_NEF_TRACEN("---------------------");
  }

};

} //namespace CGAL
#endif //CGAL_SNC_SM_OVERLAYER_H
