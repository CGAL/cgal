// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: 
// $Id: 
// 
//
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_CD3_REFLEX_EDGE_SEARCHER_H
#define CGAL_CD3_REFLEX_EDGE_SEARCHER_H

#include<CGAL/Nef_3/SNC_decorator.h>
#include<CGAL/Convex_decomposition_3/is_reflex_sedge.h>

CGAL_BEGIN_NAMESPACE

template<typename Nef_>
class Reflex_edge_searcher : public Modifier_base<typename Nef_::SNC_structure> {

  typedef Nef_                                            Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_structure          SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>              SNC_decorator;  

  typedef typename SNC_structure::Vertex_handle           Vertex_handle;
  typedef typename SNC_structure::Vertex_const_handle     Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_handle         Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle        Halffacet_handle;
  typedef typename SNC_structure::SHalfedge_handle        SHalfedge_handle;
  typedef typename SNC_structure::SHalfloop_handle        SHalfloop_handle;
  typedef typename SNC_structure::SFace_handle            SFace_handle;

  typedef typename SNC_structure::Vertex_iterator         Vertex_iterator;
  typedef typename SNC_structure::Volume_iterator         Volume_iterator;
  typedef typename SNC_structure::SHalfedge_iterator      SHalfedge_iterator;
  typedef typename SNC_structure::Shell_entry_iterator    Shell_entry_iterator;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator
                                  SHalfedge_around_svertex_circulator;

  typedef typename SNC_structure::Vector_3                Vector_3;
  typedef typename SNC_structure::Point_3                 Point_3;

  typedef typename SNC_structure::Sphere_point            Sphere_point;
  typedef typename SNC_structure::Sphere_circle           Sphere_circle;
  typedef typename SNC_structure::Sphere_segment          Sphere_segment;

  typedef typename std::deque<Halfedge_handle>           Edge_list;

 public:
  typedef typename Edge_list::iterator                    Reflex_sedge_iterator;
  typedef Edge_list                                       Container;

 public:
  Edge_list pos;
  Edge_list neg;
  Sphere_point dir;

  Reflex_edge_searcher(Sphere_point dir_in)
    : dir(dir_in) {} 

  int is_reflex_edge(Halfedge_handle e) {    
    SHalfedge_around_svertex_circulator 
      svc(e->out_sedge()), send(svc);
    int isrse = 0;
    CGAL_For_all(svc, send)
      isrse |= is_reflex_sedge(svc, dir);
    return isrse;
  }

  int is_reflex_sedge(SHalfedge_handle se) {
    return is_reflex_sedge(se, dir);
  }

  void operator()(SNC_structure& snc) {
    pos.clear();
    neg.clear();
    
    Vertex_iterator vi;
    CGAL_forall_vertices(vi, snc) {
      SHalfedge_iterator sei;
      for(sei = vi->shalfedges_begin();
	  sei != vi->shalfedges_end(); ++sei) {
	if(!sei->incident_sface()->mark()) continue;
	int isrse = CGAL::is_reflex_sedge<SNC_structure>(sei, dir);
	CGAL_NEF_TRACEN("isrse final " << sei->source()->source()->point()
			<< "->" << sei->source()->twin()->source()->point()
			<< ": " << isrse);
	if((isrse&1)==1) pos.push_back(sei->source()->twin());
	if((isrse&2)==2) neg.push_back(sei->source());
	/*
	if((isrse&2)==2) {
          sei->source()->mark()=false;
          sei->source()->twin()->mark()=false;
	}
	*/
      }
    }
  }

  /*
  void operator()(SNC_structure& snc) {
    pos.clear();
    neg.clear();    

    Reflex_edge_visitor rev(pos,neg,dir);
    SNC_decorator D(snc);
    Volume_iterator c;
    CGAL_forall_volumes(c, snc) {
      if(c->mark())
	for(Shell_entry_iterator shi=c->shells_begin(); shi!=c->shells_end(); ++shi) {
	  D.visit_shell_objects(SFace_handle(shi),rev);
	}
      }  
    }
  */
  void handle_new_edge(Halfedge_handle e) { 
    if(normalized(e->point()) == dir || 
       normalized(e->twin()->point()) == dir) {
      CGAL_error_msg( "should not happen");
      return;
    }

    if(e->twin()->source()->point() < 
       e->source()->point())
      e = e->twin();
    SHalfedge_around_svertex_circulator 
      svc(e->out_sedge()), send(svc);
    int pushed = 0;
    CGAL_For_all(svc, send) {
      int isrse = is_reflex_sedge(svc, dir);
      if(isrse == 0) continue;
      if((pushed&=1==0) && (isrse&1==1)) pos.push_back(svc->source());
      if((pushed&=2==0) && (isrse&2==2)) neg.push_back(svc->source());
      pushed |= isrse;
      if(pushed == 3)
	break;
    }
  }

  Edge_list& get_positive_redges() { return pos; }
  Edge_list& get_negative_redges() { return neg; }
  Reflex_sedge_iterator positive_redges_begin() { return pos.begin(); }
  Reflex_sedge_iterator positive_redges_end() { return pos.end(); }
  Reflex_sedge_iterator negative_redges_begin() { return neg.begin(); }
  Reflex_sedge_iterator negative_redges_end() { return neg.end(); }
};

CGAL_END_NAMESPACE
#endif // CGAL_CD3_REFLEX_EDGE_SEARCHER_H
