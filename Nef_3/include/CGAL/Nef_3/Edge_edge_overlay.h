// Copyright (c) 2007  Max-Planck-Institute Saarbruecken (Germany).
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

#ifndef CGAL_EDGE_EDGE_OVERLAY_H
#define CGAL_EDGE_EDGE_OVERLAY_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_S2/Normalizing.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_constructor.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_sphere_map.h>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 71
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <typename SNC_structure_>
class Edge_edge_overlay
{ 
public:
  typedef SNC_structure_ SNC_structure;
  typedef typename SNC_structure::Items                      Items;
  typedef typename SNC_structure_::Kernel                    Kernel;
  typedef typename Kernel::RT                                RT;
  typedef CGAL::Edge_edge_overlay<SNC_structure>             Self;
  typedef CGAL::SNC_decorator<SNC_structure>                 SNC_decorator;
  typedef typename CGAL::SNC_constructor<Items, SNC_structure>  SNC_constructor;

  typedef typename SNC_structure::Sphere_map             Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>                 SM_decorator;  
  typedef CGAL::SM_const_decorator<Sphere_map>           SM_const_decorator;
  typedef CGAL::SM_point_locator<SM_decorator>           SM_point_locator;

  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;

  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;

  typedef typename SNC_structure::Vertex_const_handle Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_handle Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle Halffacet_const_handle;
  typedef typename SNC_structure::Halffacet_const_iterator Halffacet_const_iterator;

  typedef typename SNC_structure::SVertex_iterator SVertex_iterator;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SFace_iterator SFace_iterator;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;

  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle SFace_handle;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;

  typedef typename SNC_structure::SVertex_const_handle SVertex_const_handle; 
  typedef typename SNC_structure::SHalfedge_const_handle SHalfedge_const_handle; 
  typedef typename SNC_structure::SHalfloop_const_handle SHalfloop_const_handle; 
  typedef typename SNC_structure::SFace_const_handle SFace_const_handle; 

  typedef typename SNC_structure::SHalfedge_around_facet_circulator 
    SHalfedge_around_facet_circulator;
  typedef typename SNC_structure::SHalfedge_around_facet_const_circulator 
    SHalfedge_around_facet_const_circulator;
  typedef typename SNC_structure::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename SNC_structure::SFace_cycle_const_iterator SFace_cycle_const_iterator;
  typedef typename SNC_structure::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename SNC_structure::Halffacet_cycle_const_iterator Halffacet_cycle_const_iterator;

  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Aff_transformation_3 Aff_transformation_3;

  typedef typename SNC_structure::Sphere_point Sphere_point;
  typedef typename SNC_structure::Sphere_segment Sphere_segment;
  typedef typename SNC_structure::Sphere_circle Sphere_circle;

  typedef typename SNC_structure::Mark Mark;

  typedef typename SM_decorator::SHalfedge_around_svertex_circulator 
                                 SHalfedge_around_svertex_circulator;
  typedef typename SM_decorator::SHalfedge_around_sface_circulator 
                                 SHalfedge_around_sface_circulator;
  typedef typename SM_const_decorator::SHalfedge_around_svertex_const_circulator 
                                       SHalfedge_around_svertex_const_circulator; 
  
  SNC_structure& snc;
  Halfedge_const_handle e0, e1;
  SM_decorator D;
  SM_const_decorator E;

  Edge_edge_overlay( SNC_structure& snc_,
		     Halfedge_const_handle e0_,
		     Halfedge_const_handle e1_) 
    : snc(snc_), e0(e0_), e1(e1_) {}

  template<typename Selection>
  Sphere_map* overlay_isolated_edges(const Point_3& p,
				     const Selection& BOP, bool inv) {

    Vertex_handle v = snc.new_vertex(p);
    v->mark() = BOP(e0->mark(), e1->mark(), inv);
    D = SM_decorator(&*v);
    SFace_handle sf = D.new_sface();
    sf->mark() = BOP(e0->incident_sface()->mark(),
		     e1->incident_sface()->mark(), inv);
    SVertex_handle sv = D.new_svertex(e0->point());
    sv->mark() = BOP(e0->mark(), e1->incident_sface()->mark(), inv);
    D.link_as_isolated_vertex(sv, sf);
    sv = D.new_svertex(e0->twin()->point());
    sv->mark() = BOP(e0->mark(), e1->incident_sface()->mark(), inv);
    D.link_as_isolated_vertex(sv, sf);
    sv = D.new_svertex(e1->point());
    sv->mark() = BOP(e0->incident_sface()->mark(), e1->mark(), inv);
    D.link_as_isolated_vertex(sv, sf);
    sv = D.new_svertex(e1->twin()->point());
    sv->mark() = BOP(e0->incident_sface()->mark(), e1->mark(), inv);
    D.link_as_isolated_vertex(sv, sf);
    return D.sphere_map();
  }

  template<typename Selection, typename Association>
  Sphere_map* create_edge_edge_overlay( const Point_3& p,
					const Selection& BOP, bool inv,
					Association& ) {

    //    CGAL_NEF_SETDTHREAD(43*71);
    CGAL_NEF_TRACEN(std::endl << "edge_edge " << p );

    E = SM_const_decorator(&*e1->source());
    if(E.is_isolated(e0)) {
      if(E.is_isolated(e1)) {
	return overlay_isolated_edges(p, BOP, inv);
      } else {
	std::swap(e0, e1);
	inv = !inv;
	E = SM_const_decorator(&*e1->source());
      }
    }
    
    SNC_constructor C(snc);
    Vertex_handle v(C.create_from_edge(e0,p));
    v->mark() = BOP(v->mark(), e1->mark(), inv);
    D = SM_decorator(&*v);
    SVertex_handle sv[4];
    sv[0] = sv[1] = v->svertices_begin();
    ++sv[1];

    Vector_3 vec0 = sv[0]->point() - CGAL::ORIGIN;
    Vector_3 vec1 = e1->source()->point() - CGAL::ORIGIN;
    Plane_3 mid_plane(Point_3(0,0,0), cross_product(vec0, vec1));
    Sphere_segment test_seg(sv[0]->point(),
			    CGAL::ORIGIN + mid_plane.orthogonal_vector());
    Oriented_side test_os = 
      test_seg.sphere_circle().oriented_side(e1->point());
    CGAL_assertion(test_os != ON_ORIENTED_BOUNDARY);
    if(test_os == ON_NEGATIVE_SIDE) {
      CGAL_NEF_TRACEN("change orientation of e1 " );
      e1 = e1->twin();
    }

    SM_point_locator PL(&*v);
    Object_handle o2 = PL.locate(e1->point());
    Object_handle o3 = PL.locate(e1->twin()->point());
    sv[2] = D.new_svertex(e1->point());
    sv[3] = D.new_svertex(e1->twin()->point());
    
    for(int i=0; i<4; ++i)
      CGAL_NEF_TRACEN("sv[" << i << "]= " << sv[i]->point());

    bool equator[4];
    equator[0] = equator[1] = equator[2] = equator[3] = false;
    bool on_sface[2];

    SHalfedge_handle se0, se1;
    SHalfedge_handle previous_first, previous_last;
    SFace_handle sf0, sf1;
    SHalfedge_around_svertex_circulator 
      seb[2], see[2];
    SHalfedge_around_svertex_const_circulator 
      scb[2], sce[2];
    bool empty_e[2];
    bool empty_c[2];
    empty_e[0] = empty_e[1] = empty_c[0] = empty_c[1] = false;

    on_sface[0] = CGAL::assign(sf0, o2);
    on_sface[1] = CGAL::assign(sf1, o3);

    //    set seb[0], see[1]
    if(on_sface[0]) {
      CGAL_NEF_TRACEN("found sf0 " );
      sv[2]->mark() = BOP(sf0->mark(), e1->mark(), inv);
      SFace_cycle_iterator sfci(sf0->sface_cycles_begin());
      CGAL_assertion(sfci.is_shalfedge());
      SHalfedge_handle se_tmp(sfci);
      see[1] = se_tmp;
      if(see[1]->source() != sv[0]) 
	see[1] = see[1]->snext();
      see[1] = see[1]->sprev()->twin();
      seb[0] = see[1];
    } else {
      CGAL::assign(se0, o2);
      CGAL_assertion(CGAL::assign(se0, o2));
      CGAL_NEF_TRACEN("found se0 " << se0->source()->point()
		<< "->" << se0->twin()->source()->point() );
      CGAL_NEF_TRACEN("insert sv " << sv[2]->point() );
      if(se0->source() != sv[0]) se0 = se0->twin();
      CGAL_assertion(se0->source() == sv[0]);
      sv[2]->mark() = BOP(se0->mark(), e1->mark(), inv);
      SHalfedge_handle se_new = D.split_at(se0, sv[2]);
      se_new->mark() = se_new->twin()->mark() = 
	se0->mark() = se0->twin()->mark() =
	BOP(se0->mark(), sv[2]->mark(), inv);
      // TODO: test mark of se0

      see[1] = seb[0]= se0;
      ++seb[0];
      equator[0] = equator[1] = true;
    }       

    //    set seb[1] and see[0]
    if(on_sface[1]) {
      CGAL_NEF_TRACEN("found sf1 " );
      sv[3]->mark() = BOP(sf1->mark(), e1->mark(), inv);
      SFace_cycle_iterator sfci(sf1->sface_cycles_begin());
      CGAL_assertion(sfci.is_shalfedge());
      SHalfedge_handle se_tmp(sfci);
      see[0] = se_tmp;
      if(see[0]->source() != sv[0]) 
	see[0] = see[0]->snext();
      see[0] = see[0]->sprev()->twin();
      seb[1] = see[0];
    } else {
      CGAL::assign(se1, o3);
      CGAL_assertion(CGAL::assign(se1, o3));
      CGAL_NEF_TRACEN("found se1 " << se1->source()->point()
		<< "->" << se1->twin()->source()->point() );
      CGAL_NEF_TRACEN("insert sv " << sv[3]->point() );
      if(se1->source() != sv[0]) se1 = se1->twin();
      CGAL_assertion(se1->source() == sv[0]);
      sv[3]->mark() = BOP(se1->mark(), e1->mark(), inv);
      SHalfedge_handle se_new = D.split_at(se1, sv[3]);
      se_new->mark() = se_new->twin()->mark() = 
	se1->mark() = se1->twin()->mark() =
	BOP(se1->mark(), sv[3]->mark(), inv);
      // TODO: test mark of se1

      see[0] = seb[1] = se1;
      ++seb[1];
      equator[2] = equator[3] = true;
    } 

    CGAL_assertion(seb[0]->source() == sv[0]);

    //    set empty_e if necessary
    if(seb[0] == see[0] && seb[1] == see[1]) {
      CGAL_NEF_TRACEN("both empty intervals");
      Oriented_side osx = seb[0]->circle().oriented_side(sv[3]->point());
      if(osx == ON_NEGATIVE_SIDE) {
	CGAL_NEF_TRACEN("empty 0 " );
	empty_e[0] = true;
      } else if(osx == ON_POSITIVE_SIDE) {
	empty_e[1] = true;
      } else {
	empty_e[0] = empty_e[1] = true;
      }
    } else if(seb[0] == see[0]) {
      empty_e[0] = true;
    } else if(seb[1] == see[1]) {
      empty_e[1] = true;
    }

    CGAL_NEF_TRACEN("se[01] " << (std::distance(seb[0], see[0]))
	      << ", " << (std::distance(seb[1], see[1])) );

    // finish overlay if one edge is isolated
    if(E.is_isolated(e1)) {
      Mark sfm = e1->incident_sface()->mark();
      sv[0]->mark() = BOP(sv[0]->mark(), sfm, inv);
      sv[1]->mark() = BOP(sv[1]->mark(), sfm, inv);
      SHalfedge_around_svertex_circulator 
	svc(sv[0]->out_sedge()),
	send(svc);
      CGAL_For_all(svc, send) {
	svc->incident_sface()->mark() = 
	  BOP(svc->incident_sface()->mark(), sfm, inv);
	if(svc->twin()->source() != sv[1]) continue;
	svc->mark() = BOP(svc->mark(), sfm, inv);
      }
      if(on_sface[0])
	D.link_as_isolated_vertex(sv[2], sf0);
      if(on_sface[1])
	D.link_as_isolated_vertex(sv[3], sf1);
      
      return D.sphere_map();
    }

    SHalfedge_around_svertex_const_circulator 
      svc(e1->out_sedge()), send(svc);

    // determine scb and sce for one side
    int i=0;
    bool done = false;
    Oriented_side os0 =
      svc->circle().oriented_side(sv[0]->point());
    Oriented_side os1 = os0;
    while(os0 == os1 &&
	  os1 != ON_ORIENTED_BOUNDARY &&
	  ++svc != send) 
      os1 = svc->circle().oriented_side(sv[0]->point());

    CGAL_NEF_TRACEN("osi " << os0 << ", " << os1 );
    CGAL_assertion(os1 == svc->circle().oriented_side(sv[0]->point()) ||
		   svc == send);
        
    if(os1 == ON_ORIENTED_BOUNDARY) {

      CGAL_NEF_TRACEN("svc on boundary ");

      Sphere_segment seg(sv[2]->point(), sv[3]->point(),
			 svc->circle());
      int sv_index =
	seg.has_on(sv[0]->point()) ? 0 : 1;
      equator[sv_index] = equator[sv_index+2] = true;
      
      CGAL_NEF_TRACEN("sv " << sv[0]->point() << ", " << sv[1]->point() );

      // svc is the segment from sv[2] to sv[3] that has sv[sv_index] 
      // in its interior

      if(on_sface[0]) {
	CGAL_NEF_TRACEN("se[01] " << (std::distance(seb[0], see[0]))
			<< ", " << (std::distance(seb[1], see[1])) );

	// add shalfedge_pair between sv[2] and sv[sv_index]
	SFace_cycle_iterator sfci = sf0->sface_cycles_begin();
	CGAL_assertion(sfci.is_shalfedge());
	SHalfedge_handle se_tgt(sfci);
	while(se_tgt->source() != sv[sv_index]) {
	  se_tgt = se_tgt->snext();
	  CGAL_NEF_TRACEN(se_tgt->source()->point() << " " << sv[sv_index]->point() );
	  CGAL_NEF_TRACEN(&(se_tgt->source()) << " " << &sv[sv_index] );
	}


	CGAL_NEF_TRACEN("sv[sv_index] " << sv[sv_index]->point());
	CGAL_NEF_TRACEN("after " << se_tgt->source()->point() 
		  << "->" << se_tgt->twin()->source()->point());

	SHalfedge_around_svertex_circulator se_next(se_tgt); 
	++se_next;
	CGAL_assertion(se_tgt->source() == sv[sv_index]);
	CGAL_NEF_TRACEN("new_shalfedge_pair" );
	SHalfedge_handle se_new = 
	  D.new_shalfedge_pair(sv[2], se_tgt, 1);

	se_new->mark() = se_new->twin()->mark() = 
	  BOP(sf0->mark(), svc->mark(), inv);
	se_new->circle() = normalized(Sphere_circle(se_new->source()->point(), 
						    se_new->twin()->source()->point()));
	se_new->twin()->circle() = se_new->circle().opposite();
	se_new->incident_sface() = se_new->twin()->incident_sface() = sf0;

	CGAL_NEF_TRACEN("seb[0] " << seb[0]->source()->point()
		  << "->" << seb[0]->twin()->source()->point() );
	CGAL_NEF_TRACEN("see[0] " << see[0]->source()->point()
		  << "->" << see[0]->twin()->source()->point() );

	CGAL_assertion(seb[0]->source() == sv[0]);

	CGAL_NEF_TRACEN("se[01] " << (std::distance(seb[0], see[0]))
			<< ", " << (std::distance(seb[1], see[1])) );

	if(se_next == see[0] && !empty_e[0])
	  --see[0];
	else if(se_next == see[1] && !empty_e[1])
	  --see[1];
      }

      CGAL_NEF_TRACEN("se[01] " << (std::distance(seb[0], see[0]))
		      << ", " << (std::distance(seb[1], see[1])) );
   

      if(on_sface[1]) {
	CGAL_NEF_TRACEN("sf1->mark() " << sf1->mark());

	// add shalfedge_pair between sv[3] and sv[sv_index]
	CGAL_NEF_TRACEN("seb[0] " << seb[0]->source()->point()
		  << "->" << seb[0]->twin()->source()->point() );
	CGAL_NEF_TRACEN("see[0] " << see[0]->source()->point()
		  << "->" << see[0]->twin()->source()->point() );

	SFace_cycle_iterator sfci = sf1->sface_cycles_begin();
	CGAL_assertion(sfci.is_shalfedge());
	SHalfedge_handle se_tgt(sfci);
	while(se_tgt->source() != sv[sv_index]) {
	  se_tgt = se_tgt->snext();
	  CGAL_NEF_TRACEN(se_tgt->source()->point() << " " << sv[sv_index]->point() );
	  CGAL_NEF_TRACEN(&(se_tgt->source()) << " " << &sv[sv_index] );
	}

	SHalfedge_around_svertex_circulator se_next(se_tgt); 
	++se_next;
	
	CGAL_NEF_TRACEN("sv[sv_index] " << sv[sv_index]->point());
	CGAL_NEF_TRACEN("after " << se_tgt->source()->point() 
		  << "->" << se_tgt->twin()->source()->point());

	CGAL_assertion(se_tgt->source() == sv[sv_index]);
	CGAL_NEF_TRACEN("new_shalfedge_pair" );
	SHalfedge_handle se_new =
	  D.new_shalfedge_pair(sv[3], se_tgt, 1);

	//	CGAL_assertion(se_new->snext() == se_tgt);
	CGAL_NEF_TRACEN("se_new " << se_new->source()->point() 
		  << "->" << se_new->twin()->source()->point() );
	se_new->mark() = se_new->twin()->mark() = 
	  BOP(sf1->mark(), svc->mark(), inv);
	se_new->circle() = normalized(Sphere_circle(se_new->source()->point(), 
						    se_new->twin()->source()->point()));
	se_new->twin()->circle() = se_new->circle().opposite();
	se_new->incident_sface() = 
	  se_new->twin()->incident_sface() = sf1;

	CGAL_NEF_TRACEN("seb[0] " << seb[0]->source()->point()
		  << "->" << seb[0]->twin()->source()->point() );
	CGAL_NEF_TRACEN("see[0] " << see[0]->source()->point()
		  << "->" << see[0]->twin()->source()->point() );

	CGAL_NEF_TRACEN("se[01] " << (std::distance(seb[0], see[0]))
		  << ", " << (std::distance(seb[1], see[1])) );

	if(se_next == see[0] && !empty_e[0])
	  --see[0];
	else if(se_next == see[1] && !empty_e[1])
	  --see[1];

	CGAL_NEF_TRACEN("se[01] " << (std::distance(seb[0], see[0]))
		  << ", " << (std::distance(seb[1], see[1])) );
      }

      CGAL_NEF_TRACEN("see[0] " << see[0]->source()->point()
		<< "->" << see[0]->twin()->source()->point() );
      
      if(svc == svc->twin()->snext()) {
	scb[0] = scb[1] = sce[0] = sce[1] = svc;
	empty_c[0] = empty_c[1] = true;
	done = true;
      } else {
	++svc;
	os1 = svc->circle().oriented_side(sv[0]->point());
	
	CGAL_NEF_TRACEN("++os1 " << os1 );
	
	if(os1 == ON_ORIENTED_BOUNDARY) {
	  CGAL_assertion_msg(false, "not implemented, yet"); // don't forget empty_c
	  ++svc;
	  scb[0] = scb[1] = sce[0] = sce[1] = svc;
	  if(svc != send) {
	    os1 = svc->circle().oriented_side(sv[0]->point());
	    i = os1 == ON_POSITIVE_SIDE ? 0 : 1;
	    --sce[i];
	  }
	  equator[1-sv_index] = equator[3-sv_index] = true;
	  done = true;
	} else {	
	  i = os1 == ON_POSITIVE_SIDE ? 0 : 1;
	  CGAL_NEF_TRACEN("second change found on side " << os1);
	  scb[i] = svc;
	  --svc;
	  sce[1-i] = svc;
	  --svc;
	  if(svc->circle().oriented_side(sv[0]->point()) == os1) {
	    // sedges are only on one side "
	    sce[i] = scb[1-i] = sce[1-i];
	    empty_c[1-i] = true;
	    done = true;
	  }
	}
      }
    } else if(svc == send) {
      // everything is on one side
      CGAL_NEF_TRACEN("svc == send");
      Vector_3 
	vec1(svc->source()->point() - CGAL::ORIGIN),
	vec2(svc->circle().orthogonal_vector());
      Sphere_point sp1(CGAL::ORIGIN + cross_product(vec2,vec1));
      while(svc->sprev()->twin()->circle().oriented_side(sp1) 
	    == ON_POSITIVE_SIDE) {
	++svc;
	vec1 = vec2;
	vec2 = svc->circle().orthogonal_vector();
	sp1 = CGAL::ORIGIN + cross_product(vec2,vec1);
      }
      i = os1 == ON_POSITIVE_SIDE ? 0 : 1;
      //      ++svc;
      scb[i] = sce[i] = svc;
      empty_c[1-i] = true;
      done = true;

    } else {      
      CGAL_NEF_TRACEN("found sedges on both sides ");
      //      CGAL_assertion_msg(false, "not implemented, yet");
      // Code should work, but is not tested
      CGAL_assertion(os0 != os1);
      i = os1 == ON_POSITIVE_SIDE ? 0 : 1;
      sce[1-i] = scb[i] = svc;
    }
    //    CGAL_assertion(scb[1-i] == svc);

    // determine scb, sce for other side if necessary
    if(!done) {
      CGAL_NEF_TRACEN("not done yet");
      // both sides are not empty
      os0 = svc->circle().oriented_side(sv[0]->point());
      CGAL_assertion(os0 != ON_ORIENTED_BOUNDARY);
      do {
	++svc;
	os1 = svc->circle().oriented_side(sv[0]->point());
      } while(os1 == os0);

      sce[i] = scb[1-i] = svc;
      if(os1 == ON_ORIENTED_BOUNDARY) {
	CGAL_assertion_msg(false, "degenerate case not implemented, yet");
	++scb[1-i];
	CGAL_assertion(svc->circle().has_on(sv[0]->point()) == i);
	equator[i] = equator[2+i] == true;
      }
    }

    for(int i=0; i<4; ++i)
      CGAL_NEF_TRACEN("equator[" << i << "]= " << equator[i]);
    
    CGAL_NEF_TRACEN("svmark " << sv[0]->mark() << ", " << sv[1]->mark());
    
    // TODO: marks if sv[0] and sv[1] lie on edge
    sv[0]->mark() = 
      BOP(sv[0]->mark(), sce[empty_c[0]]->twin()->incident_sface()->mark(), inv);
    sv[1]->mark() = 
      BOP(sv[1]->mark(), scb[empty_c[0]]->twin()->incident_sface()->mark(), inv);

    CGAL_NEF_TRACEN("empty[0] " << empty_c[0] << ", " << empty_e[0] );

    CGAL_assertion(seb[0]->source() == sv[0]);

    CGAL_assertion_msg(!empty_c[0] || !empty_c[1], 
		       "not implemented, yet");
    CGAL_assertion_msg(!empty_e[0] || !empty_e[1], 
		       "not implemented, yet");   

    bool first_first = true;
    bool first_last = true;
    if(!empty_c[0] && !empty_e[0]) {
      CGAL_NEF_TRACEN("sv[2] " << sv[2]->point());

      if(equator[0] || equator[1]) {
	first_first = false;
	previous_first = sv[2]->out_sedge();
	if(equator[0] && previous_first->twin()->source() != sv[0]) 
	  previous_first = previous_first->twin()->snext();
	CGAL_assertion(!equator[0] || previous_first->twin()->source() == sv[0]);
      }

      if(equator[2] || equator[3]) {
	first_last = false;
	previous_last = sv[3]->out_sedge();
	if(equator[2] && previous_last->twin()->source() != sv[0])
	  previous_last = previous_last->twin()->snext();
	CGAL_assertion(!equator[2] || previous_last->twin()->source() == sv[0]);
      }

      SHalfedge_around_svertex_const_circulator curr_outer(scb[0]);
      CGAL_For_all(curr_outer, sce[0]) {
	CGAL_NEF_TRACEN("outer " << curr_outer->incident_sface()->mark() );

	SHalfedge_handle previous_inner;
	Sphere_segment seg1(sv[2]->point(), sv[3]->point(), 
			    curr_outer->circle());
	SHalfedge_around_svertex_circulator curr_inner(seb[0]);

	CGAL_For_all(curr_inner, see[0]) {
	  CGAL_NEF_TRACEN("inner " << curr_inner->incident_sface()->mark() );
	  CGAL_assertion(!curr_inner->circle().has_on(sv[2]->point()));
	  CGAL_assertion(!curr_outer->circle().has_on(sv[0]->point()));
	  Sphere_segment seg0(sv[0]->point(), sv[1]->point(), curr_inner->circle());
	  Sphere_segment segX(curr_inner->source()->point(), 
			      curr_inner->twin()->source()->point(), 
			      curr_inner->circle());
	  CGAL_assertion(!seg0.is_long());
	  CGAL_assertion(!seg1.is_long());
	  Sphere_point sp = normalized(seg0.intersection(seg1));
	  CGAL_NEF_TRACEN("intersections oben " << sp );
	  CGAL_assertion(sp != sv[0]->point());
	  CGAL_assertion(sp != sv[1]->point());
	  CGAL_assertion(sp != sv[2]->point());
	  CGAL_assertion(sp != sv[3]->point());
	  CGAL_assertion(segX.has_on(sp));
	  CGAL_assertion(curr_inner == seb[0] || 
			 Sphere_segment(sv[2]->point(), sp).has_on
			 (previous_inner->twin()->source()->point()));
	  CGAL_assertion(seg0.source() == sv[0]->point());
	  CGAL_assertion(seg0.target() == sv[1]->point());
	  CGAL_assertion(seg1.source() == sv[2]->point());
	  CGAL_assertion(seg1.target() == sv[3]->point());
	  CGAL_assertion(seg0.has_on(sp));
	  CGAL_assertion(seg1.has_on(sp));
	  CGAL_NEF_TRACEN("seg 0 " << seg0.source() << "->" << seg0.target());
	  CGAL_NEF_TRACEN("seg 1 " << seg1.source() << "->" << seg1.target() << ":" << seg1.sphere_circle());
	  
	  SHalfedge_handle along = D.split_at(curr_inner, sp);
	  along->source()->mark() = BOP(curr_inner->mark(), 
					curr_outer->mark(), inv);
	  
	  CGAL_NEF_TRACEN("first_first " << first_first);

	  SHalfedge_handle across;
	  if(curr_inner == seb[0])
	    across = first_first ?
	      D.new_shalfedge_pair(sv[2], along, -1) :
	      D.new_shalfedge_pair(previous_first, along, -1, -1);
	  else
	    across = 
	      D.new_shalfedge_pair(previous_inner->twin(), 
				   curr_inner->twin(), -1, 1);

	  across->circle() = curr_outer->circle();
	  across->twin()->circle() = curr_outer->twin()->circle();
	  CGAL_NEF_TRACEN("across " << across->source()->point() 
		    << "->" << across->twin()->source()->point() 
		    << ":" << across->circle());
	  CGAL_NEF_TRACEN("across " << normalized(across->circle())
		    << ", " << normalized(Sphere_circle(across->source()->point(),
							across->twin()->source()->point())));
	  CGAL_assertion(normalized(across->circle()) ==
			 normalized(Sphere_circle(across->source()->point(),
						  across->twin()->source()->point())));
	  across->mark() = across->twin()->mark() = 
	    BOP(curr_inner->twin()->incident_sface()->mark(), curr_outer->mark(), inv);
	  along->mark() = along->twin()->mark() = 
	    BOP(curr_inner->mark(), curr_outer->twin()->incident_sface()->mark(), inv);
	  
	  CGAL_NEF_TRACEN("across " << across->source()->point() 
		    << "->" << across->twin()->source()->point() 
		    << ":" << across->circle() );

	  if(curr_inner == seb[0])
	    previous_first = across;
	  
	  if(!first_first) {
	    SFace_handle sf_new = D.new_sface();
	    D.link_as_face_cycle(across->twin(), sf_new);
	    sf_new->mark() = BOP(curr_inner->twin()->incident_sface()->mark(), 
				 curr_outer->twin()->incident_sface()->mark(), inv);
	    CGAL_NEF_TRACEN("new sface " << sf_new->mark());
	  }
	  previous_inner = curr_inner;
	  first_first = false;
	}
	
	previous_last = first_last ?	
	  D.new_shalfedge_pair(sv[3], previous_inner->twin(), -1) :
	  D.new_shalfedge_pair(previous_last, previous_inner->twin(), 1, -1);
	previous_last->circle() = curr_outer->twin()->circle();
	previous_last->twin()->circle() = curr_outer->circle();
	CGAL_assertion(previous_last->circle() == 
		       normalized(Sphere_circle(previous_last->source()->point(), 
						previous_last->twin()->source()->point())));
	previous_last->mark() = previous_last->twin()->mark() =
	  BOP(previous_inner->incident_sface()->mark(), curr_outer->mark(), inv);

	if(!first_last) {
	  SFace_handle sf_new = D.new_sface();
	  D.link_as_face_cycle(previous_last, sf_new);
	  sf_new->mark() = BOP(previous_inner->incident_sface()->mark(),
			       curr_outer->twin()->incident_sface()->mark(), inv);
	}
	first_last = false;
      }

      CGAL_NEF_TRACEN("sce[0]->mark() " << sce[0]->twin()->incident_sface()->mark());
      SHalfedge_around_svertex_circulator curr_inner(seb[0]);
      CGAL_For_all(curr_inner, see[0]) {
	CGAL_NEF_TRACEN("schleife oben " << curr_inner->source()->point()
			<< "->" << curr_inner->twin()->source()->point());
	curr_inner->mark() = curr_inner->twin()->mark() = 
	  BOP(curr_inner->mark(), sce[0]->twin()->incident_sface()->mark(), inv);
	if(!equator[0] && curr_inner == seb[0]) continue;
	SFace_handle sf = curr_inner->twin()->incident_sface();
	SHalfedge_around_sface_circulator hfc(curr_inner->twin()), hend(hfc);
	CGAL_For_all(hfc,hend) hfc->incident_sface() = sf;
	sf->mark() = BOP(curr_inner->twin()->incident_sface()->mark(),
			 sce[0]->twin()->incident_sface()->mark(), inv);	
      }
      CGAL_assertion(curr_inner == see[0]);

      CGAL_NEF_TRACEN("before equator[0] oben " << curr_inner->source()->point()
		<< "->" << curr_inner->twin()->source()->point());

      if(equator[2]) {
	SFace_handle sf = curr_inner->twin()->incident_sface();
	SHalfedge_around_sface_circulator hfc(curr_inner->twin()), hend(hfc);
	CGAL_For_all(hfc,hend) hfc->incident_sface() = sf;
	sf->mark() = BOP(curr_inner->twin()->incident_sface()->mark(),
			 sce[0]->twin()->incident_sface()->mark(), inv);
      } else
	--curr_inner;

      CGAL_NEF_TRACEN("final sface oben " << curr_inner->source()->point()
		<< "->" << curr_inner->twin()->source()->point());

      SFace_handle sf = curr_inner->incident_sface();
      SHalfedge_around_sface_circulator hfc(curr_inner), hend(hfc);
      CGAL_For_all(hfc,hend) hfc->incident_sface() = sf;

    } else if(empty_e[0] && !empty_c[0]) {      

      CGAL_assertion_msg(!equator[0] && !equator[1] &&
			 !equator[2] && !equator[3],
			 "not implemented, yet" );
      CGAL_assertion_msg(empty_c[1], "not implemented, yet");
      CGAL_assertion_msg(sv[2]->out_sedge() == SHalfedge_handle(),
			 "not implemented, yet");
      CGAL_assertion_msg(sv[3]->out_sedge() == SHalfedge_handle(),
			 "not implemented, yet");
      CGAL_assertion_msg(!empty_e[1], "not implemented, yet");

      Mark outer_sf = seb[1]->twin()->incident_sface()->mark();      
      SHalfedge_handle se_prev = 
	D.new_shalfedge_pair(sv[2], sv[3]);
      se_prev->circle() = scb[0]->circle();
      se_prev->twin()->circle() = scb[0]->twin()->circle();
      se_prev->mark() = se_prev->twin()->mark() = 
	BOP(outer_sf, scb[0]->mark());
      SHalfedge_around_svertex_const_circulator curr(scb[0]);
      ++curr;
      while(curr != scb[0]) {
	se_prev = D.new_shalfedge_pair(se_prev, se_prev->twin(), 1, -1);
	//	SM_io_parser<SM_decorator>::dump(D, std::cerr);
	se_prev->circle() = curr->circle();
	se_prev->twin()->circle() = curr->twin()->circle();
	se_prev->mark() = se_prev->twin()->mark() = 
	  BOP(outer_sf, curr->mark());
	SFace_handle sf_new = D.new_sface();
	sf_new->mark() = 
	  BOP(outer_sf, curr->twin()->incident_sface()->mark(), inv);
	D.link_as_face_cycle(se_prev->twin(), sf_new);
	++curr;
      }
      seb[1]->twin()->incident_sface()->mark() =
	BOP(outer_sf, scb[0]->twin()->incident_sface()->mark(), inv);
      D.link_as_face_cycle(se_prev, seb[1]->twin()->incident_sface());
      //      SM_io_parser<SM_decorator>::dump(D, std::cerr);

    } else if(empty_c[0] && !empty_e[0]) {

      CGAL_assertion_msg(!equator[0] && !equator[1] &&
			 !equator[2] && !equator[3],
			 "not implemented, yet" );
      CGAL_assertion_msg(empty_e[1], "not implemented, yet");
      CGAL_assertion_msg(!empty_c[1], "not implemented, yet");
      
      Mark outer_sf = scb[1]->twin()->incident_sface()->mark();      
      SHalfedge_around_svertex_circulator curr(seb[0]);
      CGAL_For_all(curr, see[0]) {
	curr->mark() = curr->twin()->mark() =
	  BOP(curr->mark(), outer_sf, inv);
	curr->incident_sface()->mark() =
	  BOP(curr->incident_sface()->mark(), outer_sf, inv);
      }
    }
    // else 
    //      CGAL_assertion_msg(false, "not implemented, yet");



    CGAL_NEF_TRACEN("empty[1] " << empty_c[1] << ", " << empty_e[1] );

    CGAL_NEF_TRACEN("[1] " << (std::distance(scb[1], sce[1]))
	      << ", " << (std::distance(seb[1], see[1])) );
    
    if(!empty_c[1] && !empty_e[1]) {
      CGAL_assertion(first_first == empty_c[0]);
      if(first_first) { // nothing happend on the other half
	if(equator[2] || equator[3]) {
	  first_first = false;
	  previous_first = sv[3]->out_sedge();
	  if(equator[3] && previous_first->twin()->source() != sv[1]) 
	    previous_first = previous_first->twin()->snext();
	  CGAL_assertion(!equator[3] || previous_first->twin()->source() == sv[1]);
	}

	if(equator[0] || equator[1]) {
	  first_last = false;
	  previous_last = sv[2]->out_sedge();
	  if(equator[1] && previous_last->twin()->source() != sv[1])
	    previous_last = previous_last->sprev()->twin();
	  CGAL_assertion(!equator[1] || previous_last->twin()->source() == sv[1]);
	}
      } else { // rotate over gap on unhandled half
	CGAL_assertion(!first_first && !first_last);
	std::swap(previous_last, previous_first);
	previous_first = previous_first->twin()->snext();
	if(equator[0])
	  previous_first = previous_first->twin()->snext();
	previous_last = previous_last->sprev()->twin();
	if(equator[2])
	  previous_last = previous_last->sprev()->twin();
      }

      SHalfedge_around_svertex_const_circulator curr_outer(sce[1]);
      do {
	--curr_outer;
	CGAL_NEF_TRACEN("outer " << curr_outer->incident_sface()->mark());
	CGAL_assertion(curr_outer->source()->point() == sv[2]->point());
	SHalfedge_handle previous_inner;
	Sphere_segment seg1(sv[2]->point(), sv[3]->point(), 
			    curr_outer->circle());
	SHalfedge_around_svertex_circulator curr_inner(seb[1]);
	CGAL_For_all(curr_inner, see[1]) {
	  CGAL_NEF_TRACEN("inner " << curr_inner->incident_sface()->mark());

	  CGAL_assertion(curr_inner->source() == sv[0]);
	  Sphere_segment seg0(sv[0]->point(), sv[1]->point(), curr_inner->circle());
	  Sphere_segment segX(curr_inner->source()->point(), 
			      curr_inner->twin()->source()->point(),
			      curr_inner->circle());
	  CGAL_assertion(!seg0.is_long());
	  CGAL_assertion(!segX.is_long());
	  Sphere_point sp = normalized(seg0.intersection(seg1));

	  CGAL_NEF_TRACEN("first_first " << first_first);
	  CGAL_NEF_TRACEN("first_last " << first_last);
	  CGAL_NEF_TRACEN("intersections unten " << sp );
	  CGAL_assertion(segX.has_on(sp));
	  CGAL_assertion(curr_inner == seb[1] || 
			 Sphere_segment(sv[3]->point(), sp).has_on
			 (previous_inner->twin()->source()->point()));
	  CGAL_assertion(seg0.has_on(sp));
	  CGAL_assertion(seg1.has_on(sp));
	  CGAL_NEF_TRACEN("seg 0 " << seg0.source() << "->" << seg0.target());
	  CGAL_NEF_TRACEN("seg 1 " << seg1.source() << "->" << seg1.target() << ":" << seg1.sphere_circle());
	  
	  SHalfedge_handle along = D.split_at(curr_inner, sp);
	  along->source()->mark() = BOP(curr_inner->mark(), curr_outer->mark(), inv);
	  
	  SHalfedge_handle across;
	  if(curr_inner == seb[1])
	    across = first_first ?
	      D.new_shalfedge_pair(sv[3], along, -1) :
	      D.new_shalfedge_pair(previous_first, along, 1, -1);
	  else
	    across = 
	      D.new_shalfedge_pair(previous_inner->twin(), curr_inner->twin(), -1, 1);

	  across->circle() = curr_outer->twin()->circle();
	  across->twin()->circle() = curr_outer->circle();
	  CGAL_NEF_TRACEN("across " << across->source()->point() 
		    << "->" << across->twin()->source()->point() 
		    << ":" << across->circle());
	  CGAL_NEF_TRACEN("across " << normalized(across->circle())
		    << ", " << normalized(Sphere_circle(across->source()->point(),
							across->twin()->source()->point())));
	  CGAL_assertion(normalized(across->circle()) ==
			 normalized(Sphere_circle(across->source()->point(),
						  across->twin()->source()->point())));
	  across->mark() = across->twin()->mark() = 
	    BOP(curr_inner->twin()->incident_sface()->mark(), 
		curr_outer->mark(), inv);
	  along->mark() = along->twin()->mark() = 
	    BOP(curr_inner->mark(), 
		curr_outer->incident_sface()->mark(), inv);

	  if(curr_inner == seb[1])
	    previous_first = across;

	  if(!first_first) {
	    SFace_handle sf_new = D.new_sface();
	    D.link_as_face_cycle(across->twin(), sf_new);
	    sf_new->mark() = BOP(curr_inner->twin()->incident_sface()->mark(), 
				 curr_outer->incident_sface()->mark(), inv);
	    CGAL_NEF_TRACEN("new sface" << sf_new->mark() );
	  }
	  previous_inner = curr_inner;
	  first_first = false;
	}
	
	previous_last = first_last ?	
	  D.new_shalfedge_pair(sv[2], previous_inner->twin(), -1) :
	  D.new_shalfedge_pair(previous_last, previous_inner->twin(), -1, -1);
	previous_last->circle() = curr_outer->circle();
	previous_last->twin()->circle() = curr_outer->twin()->circle();
	CGAL_assertion(normalized(previous_last->circle()) == 
		       normalized(Sphere_circle(previous_last->source()->point(), 
						previous_last->twin()->source()->point())));
	previous_last->mark() = previous_last->twin()->mark() =
	  BOP(previous_inner->incident_sface()->mark(), curr_outer->mark(), inv);


	if(!first_last) {
	  SFace_handle sf_new = D.new_sface();
	  D.link_as_face_cycle(previous_last, sf_new);
	  sf_new->mark() = BOP(previous_inner->incident_sface()->mark(),
			       curr_outer->incident_sface()->mark(), inv);
	  CGAL_NEF_TRACEN("new sface" << sf_new->mark() );
	}

	first_last = false;
      } while(curr_outer != scb[1]);

      SHalfedge_around_svertex_circulator curr_inner(seb[1]);
      CGAL_For_all(curr_inner, see[1]) {
	CGAL_NEF_TRACEN("schleife " << curr_inner->source()->point()
		  << "->" << curr_inner->twin()->source()->point());
	curr_inner->mark() = curr_inner->twin()->mark() = 
	  BOP(curr_inner->mark(), scb[1]->twin()->incident_sface()->mark(), inv);
	if(!equator[0] && curr_inner == seb[1]) continue;
	SFace_handle sf = curr_inner->twin()->incident_sface();
	SHalfedge_around_sface_circulator hfc(curr_inner->twin()), hend(hfc);
	CGAL_For_all(hfc,hend) hfc->incident_sface() = sf;
	sf->mark() = BOP(curr_inner->twin()->incident_sface()->mark(),
			 scb[1]->twin()->incident_sface()->mark(), inv);	
      }
      CGAL_assertion(curr_inner == see[1]);

      CGAL_NEF_TRACEN("before equator[2] " << curr_inner->source()->point()
		<< "->" << curr_inner->twin()->source()->point());

      if(equator[2] || 
	 (!empty_c[0] && !empty_e[0])) { 
	SFace_handle sf = curr_inner->twin()->incident_sface();
	SHalfedge_around_sface_circulator hfc(curr_inner->twin()), hend(hfc);
	CGAL_For_all(hfc,hend) hfc->incident_sface() = sf;
	sf->mark() = BOP(curr_inner->twin()->incident_sface()->mark(),
			 scb[1]->twin()->incident_sface()->mark(), inv);
      } else
	--curr_inner;

      CGAL_NEF_TRACEN("final sface " << curr_inner->source()->point()
		<< "->" << curr_inner->twin()->source()->point());

      SFace_handle sf = curr_inner->incident_sface();
      SHalfedge_around_sface_circulator hfc(curr_inner), hend(hfc);
      CGAL_For_all(hfc,hend) hfc->incident_sface() = sf;     
    } else if(empty_e[1] && !empty_c[1]) {      

      CGAL_assertion_msg(!equator[0] && !equator[1] &&
			 !equator[2] && !equator[3],
			 "not implemented, yet" );
      CGAL_assertion_msg(empty_c[0], "not implemented, yet");
      CGAL_assertion_msg(sv[2]->out_sedge() == SHalfedge_handle(),
			 "not implemented, yet");
      CGAL_assertion_msg(sv[3]->out_sedge() == SHalfedge_handle(),
			 "not implemented, yet");
      CGAL_assertion_msg(!empty_e[0], "not implemented, yet");

      Mark outer_sf = seb[0]->twin()->incident_sface()->mark();      
      SHalfedge_handle se_prev = 
	D.new_shalfedge_pair(sv[2], sv[3]);
      se_prev->circle() = scb[1]->circle();
      se_prev->twin()->circle() = scb[1]->twin()->circle();
      se_prev->mark() = se_prev->twin()->mark() = 
	BOP(outer_sf, scb[1]->mark());
      SHalfedge_around_svertex_const_circulator curr(scb[1]);
      ++curr;
      while(curr != scb[1]) {
	se_prev = D.new_shalfedge_pair(se_prev, se_prev->twin(), 1, -1);
	//	SM_io_parser<SM_decorator>::dump(D, std::cerr);
	se_prev->circle() = curr->circle();
	se_prev->twin()->circle() = curr->twin()->circle();
	se_prev->mark() = se_prev->twin()->mark() = 
	  BOP(outer_sf, curr->mark());
	SFace_handle sf_new = D.new_sface();
	sf_new->mark() = 
	  BOP(outer_sf, curr->twin()->incident_sface()->mark(), inv);
	D.link_as_face_cycle(se_prev->twin(), sf_new);
	++curr;
      }
      
      seb[0]->twin()->incident_sface()->mark() =
	BOP(outer_sf, scb[1]->twin()->incident_sface()->mark(), inv);
      D.link_as_face_cycle(se_prev, seb[0]->twin()->incident_sface());

    } else if(empty_c[1] && !empty_e[1]) {

      CGAL_assertion_msg(!equator[0] && !equator[1] &&
			 !equator[2] && !equator[3],
			 "not implemented, yet" );
      CGAL_assertion_msg(empty_e[0], "not implemented, yet");
      CGAL_assertion_msg(!empty_c[0], "not implemented, yet");
      
      Mark outer_sf = scb[0]->twin()->incident_sface()->mark();      
      SHalfedge_around_svertex_circulator curr(seb[1]);
      CGAL_For_all(curr, see[1]) {
	curr->mark() = curr->twin()->mark() =
	  BOP(curr->mark(), outer_sf, inv);
	curr->incident_sface()->mark() =
	  BOP(curr->incident_sface()->mark(), outer_sf, inv);
      }
    }
    //    else
    //      CGAL_assertion_msg(false, "not implemented, yet");

//    SM_io_parser<SM_decorator>::dump(D, std::cerr);    

    return D.sphere_map();
  }

};

} //namespace CGAL
#endif //CGAL_EDGE_EDGE_OVERLAY_H
