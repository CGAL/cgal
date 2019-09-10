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
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
//                 Micahel Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_SNC_SIMPLIFY_H
#define CGAL_SNC_SIMPLIFY_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_S2/SM_decorator.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 41
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template<typename SNC_structure>
class SNC_simplify_base : public SNC_decorator<SNC_structure> {
  
  typedef CGAL::SNC_simplify_base<SNC_structure>        Self;
  typedef CGAL::SNC_decorator<SNC_structure>            SNC_decorator;
  typedef typename SNC_structure::Sphere_map            Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>                SM_decorator;
  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Volume_handle Volume_handle;
  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;
  typedef typename SNC_structure::SFace_handle SFace_handle;

  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator; 
  typedef typename SNC_structure::Volume_iterator Volume_iterator;
  typedef typename SNC_structure::SVertex_iterator SVertex_iterator;  
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;
  typedef typename SNC_structure::SFace_iterator SFace_iterator; 

  typedef typename SNC_structure::SHalfedge_around_facet_circulator 
                                  SHalfedge_around_facet_circulator;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator 
                                  SHalfedge_around_svertex_circulator;
  typedef typename SNC_structure::SHalfedge_around_sface_circulator 
                                  SHalfedge_around_sface_circulator;

  typedef typename SNC_structure::Halffacet_cycle_iterator 
                                  Halffacet_cycle_iterator;
  typedef typename SNC_structure::SFace_cycle_iterator 
                                  SFace_cycle_iterator;

  typedef typename SNC_decorator::Sphere_kernel Sphere_kernel;
  typedef typename SNC_decorator::Sphere_point Sphere_point;
  typedef typename SNC_decorator::Sphere_segment Sphere_segment;
  typedef typename SNC_decorator::Sphere_circle Sphere_circle;
  typedef typename SNC_decorator::Sphere_direction Sphere_direction;

  typedef typename SNC_structure::Infi_box Infi_box;

  bool simplified;

 public:
  SNC_simplify_base(SNC_structure& sncs) : SNC_decorator(sncs), simplified(false) {}
  virtual ~SNC_simplify_base() {}

  typedef typename Union_find< Volume_handle>::handle UFH_volume;
  typedef typename Union_find< Halffacet_handle>::handle UFH_facet;
  typedef typename Union_find< SFace_handle>::handle UFH_sface;

  void remove_f_including_all_edge_uses_in_its_boundary_cycles
    ( Halffacet_handle f,
      Unique_hash_map< SFace_handle, UFH_sface>& hash,
      Union_find< SFace_handle>& uf )
    /* removes f and its boundary cycles, and merges up the sphere facets
       incident to them. */ {
    Halffacet_cycle_iterator fc;
    CGAL_forall_facet_cycles_of(fc, f) {
      if(fc.is_shalfedge() ) {
	SHalfedge_handle e(fc);
	SHalfedge_around_facet_circulator u(e), eend(e);

	//CGAL_For_all(u, eend) {
	for ( bool _circ_loop_flag = ! ::CGAL::is_empty_range( u, eend);
	      _circ_loop_flag;				\
	      _circ_loop_flag = ( u != eend )
	      ) {

	  SFace_handle fu = u->incident_sface(), ftu = u->twin()->incident_sface();
	  //	  CGAL_NEF_TRACEN("sfUNION of "<<IO->index(fu)<<" & "<<IO->index(ftu));
	  merge_sets( fu, ftu, hash, uf);
	  SM_decorator SD(&*u->source()->source());
	  //	  CGAL_NEF_TRACEN("removing "<<IO->index(u)<<" & "<<IO->index(u->twin()));
	  Halfedge_handle src(u->source()), tgt(u->target());
	  if ( SD.is_closed_at_source(u) ) 
	    SD.set_face( src, fu);
	  if ( SD.is_closed_at_source( u->twin()) ) 
 	    SD.set_face( tgt, fu);
	  /* TO VERIFY: does is_closed_at_source(u) imply is_isolated(src)?
	     if it is true, the svertex face update is not necesary. */

	  SHalfedge_around_facet_circulator next = u;
	  ++next;
	  SD.delete_edge_pair(u); 
	  if( SD.is_isolated(src))
	    SD.set_face(src,fu);
	  SM_decorator SD2(&*tgt->source());
	  if( SD2.is_isolated(tgt))
	    SD2.set_face(tgt,fu);
	  u = next;
	  /* TO VERIFY: can both svertices be isolated at the same time? */
   }
      }
      else if(fc.is_shalfloop()) {
	SHalfloop_handle l(fc);
	// this code is currenlty not used, but it is potentially need 
	// in the future, e.g for complex marks or a relative interior 
	// function 
	SFace_handle fu = l->incident_sface(), ftu = l->twin()->incident_sface();
	//	CGAL_NEF_TRACEN("UNION of "<<IO->index(fu)<<" & "<<IO->index(ftu));
	merge_sets( fu, ftu, hash, uf);
	SM_decorator SD(&*l->incident_sface()->center_vertex());
	//	CGAL_NEF_TRACEN("removing "<<IO->index(l)<<" & "<<IO->index(l->twin()));
	SD.delete_loop_only();
      }
    }
    //    CGAL_NEF_TRACEN("removing "<<IO->index(f)<<" & "<<IO->index(f->twin()));
    this->sncp()->delete_halffacet_pair(f);
    return;
  }

  bool is_part_of_volume(Vertex_handle v) 
    /* determines if a vertex v is part of a volume, cheking if its local
       graph is trivial (only one sface with no boundary). */  {
    SM_decorator SD(&*v);
    CGAL_assertion( !is_empty_range( SD.sfaces_begin(), SD.sfaces_end()));
    if( is_empty_range( SD.svertices_begin(), SD.svertices_end()) &&
	is_empty_range( SD.shalfedges_begin(), SD.shalfedges_end()) &&
	!SD.has_shalfloop())
      return true;
    return false;
  }

  bool is_part_of_facet(Vertex_handle v) 
    /* determines if a vertex v is part of a the relative interior of a 
       facet, checking if its local graph consists just of a sloop and
       two incident sfaces. */ {
    SM_decorator SD(&*v);
    CGAL_assertion( !is_empty_range( SD.svertices_begin(),
					  SD.svertices_end()) ||
			 is_empty_range( SD.shalfedges_begin(),
					 SD.shalfedges_end()));
    return( SD.has_shalfloop() &&
	    is_empty_range( SD.svertices_begin(), SD.svertices_end()));
  }

  bool is_part_of_edge(Vertex_handle v) {
    /* determines if a vertex v is part of a edge, checking at its local 
       graph for exactly two antipodal vertices  */

    SM_decorator SD(&*v);
    if(SD.has_shalfloop())
      return false;
    if(SD.svertices_begin() == SD.svertices_end())
      return false;
    if(++(SD.svertices_begin()) == SD.svertices_end())
      return false;
    
    CGAL_NEF_TRACE(v->point()<<" is in edge interior? ");
    SVertex_iterator sv(SD.svertices_begin());
    SVertex_handle p1(sv++), p2(sv++);
    if( sv != SD.svertices_end())
      return false;
    
    CGAL_NEF_TRACE("has two svertices ");
    Sphere_point sp1(p1->point()), sp2(p2->point());
    return (sp1 == sp2.antipode());
  }

  bool simplify_redundant_box_vertex(Vertex_handle v, bool snc_computed) {
    CGAL_warning("altered code");
    return false;
    if(snc_computed) return false;
    if(!Infi_box::is_redundant_box_vertex(*v)) return false;
    this->sncp()->delete_vertex(v);
    simplified = true;
    return true;
  }

  bool simplify_redundant_vertex_in_volume(Vertex_handle v) {
    if(is_part_of_volume(v)) {
      //      CGAL_NEF_TRACEN("mark("<<IO->index(v)<<")="<<v->mark()<<", "<<
      //		      "mark("<<IO->index(v->sfaces_begin()->volume())<<")="<<
      //		      v->sfaces_begin()->mark());
      if(v->mark() == v->sfaces_begin()->mark()) {
	//	CGAL_NEF_TRACEN("removing isolated vertex "<<IO->index(v));
	this->sncp()->delete_vertex(v);
	simplified = true;
      }
      return true;
    }
    return false;
  }

  bool simplify_redundant_vertex_on_facet(Vertex_handle v) {
    if( is_part_of_facet(v)) {
      if( v->mark() == v->shalfloop()->mark()) {
	//	CGAL_NEF_TRACEN("removing "<<IO->index(v)<<
	//			", " << v->point() << " on facet ");
	this->sncp()->delete_vertex(v);
	simplified = true;
      }
      return true;
    }
    return false;
  }

  bool simplify_redundant_vertex_on_edge(Vertex_handle v, bool snc_computed) {
    if( is_part_of_edge(v)) {
      SVertex_iterator sv(v->svertices_begin());
      Halfedge_handle e1(sv++), e2(sv++);
      CGAL_assertion( sv == v->svertices_end());
      if( e1->mark() == v->mark() && e1->mark() == e2->mark()) {
	//	CGAL_NEF_TRACEN("merging "<<IO->index(e1)<<" & "<<IO->index(e2)<<
	//			" in "<<IO->index(v));
	if(snc_computed)
	  merge_halfedge_pairs( e1, e2);
	else
	    this->sncp()->delete_vertex(v);
	simplified = true;
      }
      return true;
    }
    return false;
  }

  bool simplify_redundant_vertex_on_edgeI(Vertex_handle v) {
    if( is_part_of_edge(v)) {
      SVertex_iterator sv(v->svertices_begin());
      Halfedge_handle e1(sv++), e2(sv++);
      CGAL_assertion( sv == v->svertices_end());
      if( e1->mark() == v->mark() && e1->mark() == e2->mark()) {
	SHalfedge_around_svertex_circulator 
	  sec(e1->out_sedge()), send(sec);
	CGAL_For_all(sec,send) {
	  sec->next()->prev() = sec->prev();
	  sec->prev()->next() = sec->next();
	  sec->twin()->next()->prev() = sec->twin()->prev();
	  sec->twin()->prev()->next() = sec->twin()->next();
	}
	if(e2->get_index() < e1->get_index())
	  e1->twin()->set_index(e2->twin()->get_index());
	else
	  e2->twin()->set_index(e1->twin()->get_index());	
	e1->twin()->twin() = e2->twin();
	e2->twin()->twin() = e1->twin();
	this->sncp()->delete_vertex(v);
      }
      return true;
    }
    return false;
  }

  bool vertex_in_volume_simplification() {
    simplified = false;

    Vertex_iterator v = (*this->sncp()).vertices_begin();
    while( v != (*this->sncp()).vertices_end()) {
      Vertex_iterator v_next(v);
      ++v_next;
      simplify_redundant_vertex_in_volume(v);
      v = v_next;
    }
    return simplified;
  }

  bool vertex_simplification(bool snc_computed = true) {
    simplified = false;

    Vertex_iterator v = (*this->sncp()).vertices_begin();
    while( v != (*this->sncp()).vertices_end()) {
      Vertex_iterator v_next(v);
      ++v_next;
      if(!simplify_redundant_box_vertex(v, snc_computed))
	if(!simplify_redundant_vertex_in_volume(v))
	  if(!simplify_redundant_vertex_on_facet(v))
	    simplify_redundant_vertex_on_edge(v, snc_computed);
      v = v_next;
    }
    return simplified;
  }

  void vertex_simplificationI() {
    simplified = false;

    Vertex_iterator v = (*this->sncp()).vertices_begin();
    while( v != (*this->sncp()).vertices_end()) {
      Vertex_iterator v_next(v);
      ++v_next;
      if(!simplify_redundant_vertex_in_volume(v))
	if(!simplify_redundant_vertex_on_facet(v))
	  simplify_redundant_vertex_on_edgeI(v);
      v = v_next;
    }
  }
  
  bool simplify() {

    bool update_facets  =  false;
    bool update_sfaces  =  false;
    bool update_volumes =  false;

    CGAL_NEF_TRACEN(">>> simplifying");
    SNC_decorator D(*this->sncp());
    
    Unique_hash_map< Volume_handle, UFH_volume> hash_volume;
    Unique_hash_map< Halffacet_handle, UFH_facet> hash_facet;
    Unique_hash_map< SFace_handle, UFH_sface> hash_sface;
    Union_find< Volume_handle> uf_volume;
    Union_find< Halffacet_handle> uf_facet;
    Union_find< SFace_handle> uf_sface;

    /* We discard  the information about boundary entry points, first
       on volumes, facets on sfacets.  Since during the volumes simplification
       is required the remotion of facet cycles, the information about those
       cycles is keep until the this simplification step is performed. */

    this->sncp()->clear_boundary();

    Volume_iterator c;
    CGAL_forall_volumes( c, *this->sncp()) {
      hash_volume[c] = uf_volume.make_set(c);
      this->sncp()->reset_object_list(c->shell_entry_objects());
    }
    SFace_iterator sf;
    CGAL_forall_sfaces( sf, *this->sncp()) {
      hash_sface[sf] = uf_sface.make_set(sf);
      this->sncp()->reset_sm_object_list(sf->boundary_entry_objects());
    }

    /* 
     * Volumes simplification 
     */

    Halffacet_handle f(D.halffacets_begin());
    while( f != D.halffacets_end() && f->is_twin())
      f++;
    while( f != D.halffacets_end()) {
      CGAL_assertion( !f->is_twin());
      Halffacet_iterator f_next(f);
      do
	++f_next;
      while( f_next != D.halffacets_end() && f_next->is_twin());
      CGAL_assertion( f != f->twin());
      Volume_handle c1 = f->incident_volume(), c2 = f->twin()->incident_volume();
      CGAL_NEF_TRACEN(" mark(c1)="<<c1->mark()<<
           	      " mark(f)="<<f->mark() <<
		      " mark(c2)="<<c2->mark()<<
      	              " is_f->twin()="<<f->is_twin());
      if( c1->mark() == f->mark() && c1->mark() == c2->mark()
	  && D.is_standard(f)) {
	merge_sets( c1, c2, hash_volume, uf_volume);
	remove_f_including_all_edge_uses_in_its_boundary_cycles
	  (f, hash_sface, uf_sface);
	update_sfaces = update_volumes = true;
	CGAL_NEF_TRACEN("UNION of c1 & c2");
      }
      f = f_next;
    }
    
    CGAL_forall_halffacets( f, *this->sncp()) {
      hash_facet[f] = uf_facet.make_set(f);
      this->sncp()->reset_object_list(f->boundary_entry_objects());
    }
    
    /* 
     * Edges simplification
     */

    Halfedge_iterator e(D.halfedges_begin());
    while( e != D.halfedges_end() && e->is_twin())
      e++;
    while( e != (*this->sncp()).halfedges_end()) {
      CGAL_assertion( !e->is_twin());
      Halfedge_iterator e_next(e);
      do 
	e_next++;
      while( e_next != D.halfedges_end() && e_next->is_twin());
      SM_decorator SD(&*e->source());
      if( SD.is_isolated(e)) {
	if(e->mark() == e->incident_sface()->volume()->mark()) {
	  CGAL_NEF_TRACEN("removing pair ");
	  this->sncp()->delete_halfedge_pair(e);
	  update_facets = true;
	}
      } 
      else { 
	if( D.has_outdeg_two(e)) {
	  SHalfedge_handle e1(SD.first_out_edge(e)); 
	  SHalfedge_handle e2(SD.cyclic_adj_succ(e1));
	  if( e1->circle()==e2->twin()->circle() &&
	      e1->mark()==e->mark() && e1->mark()==e2->mark()) {
	    Halffacet_handle f1(e1->facet()); 
	    Halffacet_handle f2(e2->facet());
	    CGAL_NEF_TRACEN("UNION of f1 & f2->twin()");
	    merge_sets( f1, f2->twin(), hash_facet, uf_facet);
	    merge_sets( f1->twin(), f2, hash_facet, uf_facet);
	    CGAL_NEF_TRACEN("removing e");
	    remove_edge_and_merge_facet_cycles(e);
	    update_facets = true;
	  }
	}
      }
      e = e_next;
    }

    update_facets = vertex_simplification() || update_facets;
    purge_no_find_objects(hash_volume, hash_facet, hash_sface, uf_volume, 
			  uf_facet, uf_sface);
    create_boundary_links_forall_sfaces( hash_sface, uf_sface);
    create_boundary_links_forall_facets( hash_facet, uf_facet);
    create_boundary_links_forall_volumes( hash_volume, uf_volume);
    
    CGAL_NEF_TRACEN(">>> simplifying done ");

    return update_sfaces || update_facets || update_volumes;
  }
   
  void remove_edge_and_merge_facet_cycles( Halfedge_handle e) {
    CGAL_assertion_code(SNC_decorator D(*this->sncp()));
    CGAL_assertion( D.has_outdeg_two(e));
    Halfedge_handle et = e->twin();
    CGAL_assertion( D.has_outdeg_two(et));
    CGAL_NEF_TRACEN("source " << e->source()->point());
    CGAL_NEF_TRACEN("target " << et->source()->point());
    SHalfedge_handle e1 = e->out_sedge();
    SHalfedge_handle e2 = e1->prev()->snext();
    merge_sedges_at_target_and_remove_svertex( e1->twin(), e);
    merge_sedges_at_target_and_remove_svertex( e2->twin(), et);
  }

   void merge_sedges_at_target_and_remove_svertex( SHalfedge_handle s1,
						   SVertex_handle v) {
     SNC_decorator D(*this->sncp());
     SM_decorator SD(&*v->source());
     CGAL_assertion( s1->twin()->source() == v);
     SHalfedge_handle s2(s1->snext());
     CGAL_assertion( s2->source() == v);
     //     CGAL_NEF_TRACEN("s1 = " << IO->index(s1));
     //     CGAL_NEF_TRACEN("s2 = " << IO->index(s2));
     if( s1 == s2) {
       //       CGAL_NEF_TRACEN(IO->index(s1)<<'('<<IO->index(s2->twin())<<") to sloop");
       Halffacet_handle f = s1->facet();
       // Halffacet_handle facet = s1->facet();
       SD.convert_edge_to_loop(s1);
       CGAL_assertion(SD.shalfloop() != SHalfloop_handle());
       D.add_sloop_to_facet( SD.shalfloop(), f);
       //       CGAL_NEF_TRACEN(IO->index(s2)<<" removed");
     }
     else {
       CGAL_assertion( D.has_outdeg_two(v));
       D.link_as_prev_next_pair( s1, s2->next());
       //       CGAL_NEF_TRACEN(IO->index(s1)<<" "<<IO->index(s2->next())<<" linked.");
       D.link_as_prev_next_pair(s2->next()->twin(), s1->twin());
       //       CGAL_NEF_TRACEN(IO->index(s2->next()->twin())<<" "<<
       //	      IO->index(s1->twin())<<" linked.");
       SD.merge_edge_pairs_at_target( s1); // s2 is removed
       //       CGAL_NEF_TRACEN(IO->index(s2)<<" removed");
     }
   }

   virtual void merge_halfedge_pairs( SVertex_handle p, SVertex_handle q) {
     SNC_decorator D(*this->sncp());
     CGAL_assertion( p->source() == q->source());
     Vertex_handle v(p->source()); 
     CGAL_assertion( is_part_of_edge(v));
     SM_decorator SD(&*v);
     SHalfedge_around_svertex_circulator s(SD.first_out_edge(p)), se(s);
     CGAL_For_all( s, se) {
       D.link_as_prev_next_pair( s->prev(), s->next());
       D.link_as_prev_next_pair( s->twin()->prev(), s->twin()->next());
     }
     D.make_twins( p->twin(), q->twin());
     SD.delete_vertex(p);
     SD.delete_vertex(q);
     this->sncp()->delete_vertex(v);
   }

   void purge_no_find_objects( 
      Unique_hash_map< Volume_handle, UFH_volume>& hash_volume,
      Unique_hash_map< Halffacet_handle, UFH_facet>& hash_facet,
      Unique_hash_map< SFace_handle, UFH_sface>& hash_sface,
      Union_find< Volume_handle>& uf_volume,
      Union_find< Halffacet_handle>& uf_facet,
      Union_find< SFace_handle>& uf_sface ) {
     SNC_decorator D(*this->sncp());
     SFace_iterator sf;
     std::list<SFace_handle> sflist;
     CGAL_forall_sfaces( sf, *this->sncp()) {
       if( uf_sface.find(hash_sface[sf]) != hash_sface[sf]) {
	 //	 CGAL_NEF_TRACEN("no find object "<<IO->index(sf));
	 sflist.push_back(sf);
       }
     }

     typename std::list<SFace_handle>::const_iterator sfli;
     for(sfli = sflist.begin(); sfli != sflist.end(); sfli++){     
       SM_decorator SD(&*(*sfli)->center_vertex());
       SD.delete_face_only(*sfli);
     }

     Halffacet_iterator f;
     std::list<Halffacet_handle> flist;
     CGAL_forall_facets( f, *this->sncp()) {
       //       CGAL_NEF_TRACEN("facet "<<IO->index(f));
       if( uf_facet.find(hash_facet[f]) != hash_facet[f]) {
	 //	 CGAL_NEF_TRACEN("no find object "<<IO->index(f));
	 flist.push_back(f);
       }
     }
     
     typename std::list<Halffacet_handle>::const_iterator fli;
     for(fli = flist.begin(); fli != flist.end(); fli++)
       this->sncp()->delete_halffacet_pair(*fli);

     Volume_iterator c;
     std::list<Volume_handle> clist;
     CGAL_forall_volumes( c, *this->sncp()) {
       if( uf_volume.find(hash_volume[c]) != hash_volume[c]) {
	 //	 CGAL_NEF_TRACEN("no find object "<<IO->index(c));
	 clist.push_back(c);
       }
     }

     typename std::list<Volume_handle>::const_iterator cli;
     for(cli = clist.begin(); cli != clist.end(); cli++){     
       this->sncp()->delete_volume(*cli);
     }
   }

  void create_boundary_links_forall_sfaces(
      Unique_hash_map< SFace_handle, UFH_sface>& hash,
      Union_find< SFace_handle>& uf ) {
    Unique_hash_map< SHalfedge_handle, bool> linked(false);
    SNC_decorator D(*this->sncp());
    SHalfedge_iterator e;
    CGAL_forall_shalfedges(e, *this->sncp()) {
      if( linked[e])
	continue;
      SM_decorator SD(&*e->source()->source());
      SFace_handle sf = *(uf.find(hash[e->incident_sface()]));
      CGAL_assertion( sf != SFace_handle());
      SHalfedge_around_sface_circulator c(e), cend(c);
      CGAL_For_all( c, cend) {
	SD.set_face(c, sf);
	linked[c] = true;
      }
      SD.store_sm_boundary_object( e, sf);
    }

    SVertex_handle sv;
    CGAL_forall_svertices(sv, *this->sncp()) {
      SM_decorator SD(&*sv->source());
      if( SD.is_isolated(sv)) {
	SFace_handle sf = *(uf.find(hash[sv->incident_sface()])); 
	CGAL_assertion( sf != SFace_handle());
	SD.set_face( sv, sf);
	SD.store_sm_boundary_object( sv, sf);
      }
    }

    SHalfloop_handle sl;
    CGAL_forall_shalfloops(sl, *this->sncp()) {
      SM_decorator SD(&*sl->incident_sface()->center_vertex());
      //I added the following 'if' because even if the map has been cleared, when merging edges
      //if one sloop is created and the vertex not simplified, then the map is updated.
      //In the example below, the edge e is removed and a sloop is create on the sphere map
      //of a vertex v. But since three edges are still incident to v, v is not simplified. 
      //There is one point where the distance between two non-edge-adjacent facets is 0.
      //
      //      |               |
      // |    |  \  |  / |    |
      // |    /   \ | /  |    /
      // |   /     \|/   |   /
      // |  /      /     |  /
      // | /      / e    | /
      // |/______/_______|/
      //   
      if (!SD.map()->is_sm_boundary_object(sl)) 
        SD.store_sm_boundary_object( sl, sl->incident_sface());
    }
  }

  void create_boundary_links_forall_facets(
      Unique_hash_map< Halffacet_handle, UFH_facet>& hash,
      Union_find< Halffacet_handle>& uf) {
    Unique_hash_map< SHalfedge_handle, bool> linked(false);
    SNC_decorator D(*this->sncp());
    SHalfedge_iterator u;
    CGAL_forall_shalfedges(u, *this->sncp()) {
      if( linked[u])
	continue;
      /* set find(f) as incident facet of every edge use on the cycle of u */
      SHalfedge_handle u_min = u;
      Halffacet_handle f = *(uf.find(hash[u->facet()]));
      SHalfedge_around_facet_circulator c(u), cend(c);
      CGAL_For_all( c, cend) {
	D.set_facet( c, f);
	if( lexicographically_xyz_smaller(c->source()->source()->point(), 
					  u_min->source()->source()->point()))
	  u_min = c;
	linked[c] = true;
      }
      /* store the edge use at the lexicographicaly minimum facet vertex, as
	 a cycle entry of f.  The outermost cycle is stored at first
	 on the facet's cycles list. */
      if( is_empty_range( f->boundary_entry_objects().begin(),
			  f->boundary_entry_objects().end())) {
	D.store_boundary_object( u_min, f);
	CGAL_NEF_TRACEN("new outer cycle min. vertex: "<< u_min->source()->source()->point());
      }
      else {
	SHalfedge_handle f_sedge;
	CGAL_assertion( CGAL::assign( f_sedge, 
				     f->boundary_entry_objects().front()));
	CGAL::assign( f_sedge, f->boundary_entry_objects().front());
	if( lexicographically_xyz_smaller(u_min->source()->source()->point(), 
					  f_sedge->source()->source()->point()))
	  D.store_as_first_boundary_object( u_min, f);
	else
	  D.store_boundary_object( u_min, f);
      }
    }
    SHalfloop_iterator l;
    CGAL_forall_shalfloops( l, *this->sncp()) {
      Halffacet_handle f = *(uf.find(hash[l->facet()]));
      D.set_facet( l, f);
      D.store_boundary_object( l, f);
    }
  }

  void create_boundary_links_forall_volumes( 
      Unique_hash_map< Volume_handle, UFH_volume>& hash,
      Union_find< Volume_handle>& uf) {
    typedef typename SNC_decorator::template Shell_volume_setter<SNC_decorator> Volume_setter;
    //   typedef Unique_hash_map< SFace_handle, bool> SFace_map;
    //  SFace_map linked(false);

    SNC_decorator D(*this->sncp());
    Volume_setter setter(D);

    SFace_iterator sf;
    Volume_handle c;
    CGAL_forall_sfaces(sf, *this->sncp()) {
      //      CGAL_NEF_TRACEN("SFace " << IO->index(sf));
      if( setter.is_linked(sf)) continue;
      c = *(uf.find(hash[sf->volume()]));
      //      CGAL_NEF_TRACEN("Volume " << IO->index(c));
      setter.set_volume(c);
      D.visit_shell_objects( sf, setter );      
      D.store_boundary_object( sf, c);
    }
  }
};


template<typename Items, typename SNC_structure>
class SNC_simplify : public SNC_simplify_base<SNC_structure> {
public:
  SNC_simplify(SNC_structure& sncs)
    : SNC_simplify_base<SNC_structure>(sncs)
  {}
};


template<typename SNC_structure>
class SNC_simplify<SNC_indexed_items, SNC_structure> 
 : public SNC_simplify_base<SNC_structure> {

  typedef SNC_simplify_base<SNC_structure>              Base;
  typedef CGAL::SNC_decorator<SNC_structure>            SNC_decorator;
  typedef typename SNC_structure::Sphere_map            Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>                SM_decorator;
  
  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator; 
  typedef typename SNC_structure::Halffacet_cycle_iterator 
                                  Halffacet_cycle_iterator;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator 
                                  SHalfedge_around_svertex_circulator;
  typedef typename SNC_structure::SHalfedge_around_facet_circulator 
                                  SHalfedge_around_facet_circulator;

  using Base::is_part_of_edge;

 public:
  SNC_simplify(SNC_structure& sncs) : Base(sncs) {}

  bool simplify() {

    bool result = Base::simplify();

    Halffacet_iterator fit;
    CGAL_forall_halffacets(fit, *this) {
      Halffacet_cycle_iterator fci = fit->facet_cycles_begin();
      SHalfedge_handle se_first = fci;
      int index = se_first->get_index();
      for(; fci!=fit->facet_cycles_end(); ++fci) {
	if(fci.is_shalfedge()) {
	  SHalfedge_around_facet_circulator sc1(fci), sc2(sc1);
	  CGAL_For_all(sc1,sc2)
	    sc1->set_index(index);
	} else if(fci.is_shalfloop()) {
	  SHalfloop_handle sl(fci);
	  sl->set_index(index);
	} else
	  CGAL_error_msg( "wrong handle");
      }
    }
    
    return result;
  }

  void merge_halfedge_pairs( SVertex_handle p, SVertex_handle q) {
    SNC_decorator D(*this->sncp());
    CGAL_assertion( p->source() == q->source());
    Vertex_handle v(p->source()); 
    CGAL_assertion( is_part_of_edge(v));
    SM_decorator SD(&*v);
    SHalfedge_around_svertex_circulator s(SD.first_out_edge(p)), se(s);
    CGAL_For_all( s, se) {
      D.link_as_prev_next_pair( s->prev(), s->next());
      D.link_as_prev_next_pair( s->twin()->prev(), s->twin()->next());
    }
    D.make_twins( p->twin(), q->twin());
    if(q->get_index() < p->get_index())
      p->twin()->set_index(q->twin()->get_index());
    else
      q->twin()->set_index(p->twin()->get_index());
    SD.delete_vertex(p);
    SD.delete_vertex(q);
    this->sncp()->delete_vertex(v);
  }
};

} //namespace CGAL
#endif // CGAL_SNC_STRUCTURE_H
