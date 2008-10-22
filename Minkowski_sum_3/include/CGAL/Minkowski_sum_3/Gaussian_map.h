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
#ifndef CGAL_MS3_GAUSSIAN_MAP
#define CGAL_MS3_GAUSSIAN_MAP

#include <CGAL/Nef_S2/SM_items.h>
#include <CGAL/Nef_S2/Sphere_map.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_S2/SM_io_parser.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/Minkowski_sum_3/PointMark.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 223
#include <CGAL/Nef_2/debug.h>

CGAL_BEGIN_NAMESPACE


template <class K, class Nef, class Mark_ = PointMark<K> >
class Gaussian_map :
  public CGAL::SM_decorator<CGAL::Sphere_map<CGAL::Sphere_geometry<K>, 
					     CGAL::SM_items, Mark_> > {

  typedef K                                               Kernel;
  typedef CGAL::Sphere_geometry<K>                        Sphere_kernel;
  typedef typename Kernel::Point_3                        Point_3;
  typedef typename Kernel::Vector_3                       Vector_3;
  typedef Mark_                                           Mark;
  typedef CGAL::Sphere_map<Sphere_kernel, 
                           CGAL::SM_items, Mark>          Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>                  SM_decorator;
  typedef SM_decorator                                    Base;
  typedef CGAL::SM_overlayer<SM_decorator>                SM_overlayer;
  typedef typename Sphere_kernel::Sphere_circle           Sphere_circle;
 public:
  typedef typename Sphere_map::SVertex_handle             SVertex_handle;
  typedef typename Sphere_map::SHalfedge_handle           SHalfedge_handle;
  typedef typename Sphere_map::SHalfloop_handle           SHalfloop_handle;
  typedef typename Sphere_map::SFace_handle               SFace_handle;
  typedef typename Sphere_map::SFace_iterator             SFace_iterator;
  typedef typename Sphere_map::SHalfedge_iterator         SHalfedge_iterator;
  typedef typename Sphere_map::SVertex_iterator           SVertex_iterator;
  typedef typename Sphere_map::SVertex_const_iterator     SVertex_const_iterator;
  typedef typename Sphere_map::SVertex_const_handle       SVertex_const_handle;
  typedef typename Sphere_map::SFace_cycle_iterator       SFace_cycle_iterator;
  typedef typename Sphere_map::SHalfedge_around_svertex_circulator
    SHalfedge_around_svertex_circulator;
  typedef typename Sphere_map::SHalfedge_around_sface_circulator
    SHalfedge_around_sface_circulator;
  typedef typename Sphere_map::SHalfedge_around_sface_const_circulator
    SHalfedge_around_sface_const_circulator;

  typedef typename Sphere_map::Object_handle              Object_handle;

  template<typename Nef_polyhedron_3>
  class SVertex_creator2 {

    typedef typename Nef_polyhedron_3::Vertex_const_handle      Vertex_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle   Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle    Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle   SHalfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfloop_const_handle   SHalfloop_const_handle;
    typedef typename Nef_polyhedron_3::SFace_const_handle       SFace_const_handle;

    typedef CGAL::Unique_hash_map<Halffacet_const_handle, SVertex_handle> Facet2SVertex_hash;

    SM_decorator SM;
    Facet2SVertex_hash& Facet2SVertex;
    
  public:
    SVertex_creator2(Sphere_map* smap, Facet2SVertex_hash& F2SV)
      : SM(smap), Facet2SVertex(F2SV) {}

    void visit(Vertex_const_handle v) {}
    void visit(Halfedge_const_handle e) {}
    void visit(SHalfedge_const_handle se) {}
    void visit(SHalfloop_const_handle sl) {}
    void visit(SFace_const_handle sf) {}
    void visit(Halffacet_const_handle f) {

      CGAL_NEF_TRACEN( "SVertex_creator2 " << f->twin()->plane() );
      SVertex_handle sv
	(SM.new_svertex(normalized(f->twin()->plane().orthogonal_vector())));
      sv->mark() = Mark(Point_3(0,0,0), f->mark());
      Facet2SVertex[f] = sv;
    }

  };

  template<typename Nef_polyhedron_3>
  class SVertex_creator {     

    typedef typename Nef_polyhedron_3::Vertex_const_handle      Vertex_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle   Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle    Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle   SHalfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfloop_const_handle   SHalfloop_const_handle;
    typedef typename Nef_polyhedron_3::SFace_const_handle       SFace_const_handle;

    typedef typename Nef_polyhedron_3::Halffacet_cycle_const_iterator   
      Halffacet_cycle_const_iterator;
    typedef typename Nef_polyhedron_3::SHalfedge_around_facet_const_circulator
      SHalfedge_around_facet_const_circulator;

    typedef typename Nef_polyhedron_3::Point_3                  Point_3;
    typedef typename Nef_polyhedron_3::Sphere_point             Sphere_point;

    typedef CGAL::Unique_hash_map<Halffacet_const_handle, SVertex_handle> Facet2SVertex_hash;
    typedef CGAL::Unique_hash_map<Halffacet_const_handle, bool> Facet2bool_hash;
    typedef CGAL::Unique_hash_map<Vertex_const_handle, bool> Vertex2bool_hash;
    typedef CGAL::Unique_hash_map<Halfedge_const_handle, bool> Edge2bool_hash;
    typedef CGAL::Unique_hash_map<SHalfedge_const_handle, SHalfedge_const_handle> SEdge2SEdge_hash;
    typedef CGAL::Unique_hash_map<SFace_const_handle, bool> SFace2bool_hash;

    SM_decorator SM;
    Facet2SVertex_hash& Facet2SVertex;
    Facet2bool_hash& omit_facet;
    SEdge2SEdge_hash& next;
    Vertex2bool_hash& omit_vertex;
    Edge2bool_hash& omit_edge;
    SFace2bool_hash& Shell;

  public:
    SVertex_creator(Sphere_map* smap, Facet2SVertex_hash& F2SV, Facet2bool_hash& F2b, SEdge2SEdge_hash& SE2SE, Vertex2bool_hash& V2b, Edge2bool_hash& E2b, SFace2bool_hash& SHELL)
      : SM(smap), Facet2SVertex(F2SV), omit_facet(F2b), next(SE2SE), omit_vertex(V2b), omit_edge(E2b), Shell(SHELL) {}

  private:
    bool svertex_exists(Sphere_point sp, SVertex_handle& sv) {
      SVertex_iterator svi;
      sp = normalized(sp);
      CGAL_forall_svertices(svi, SM) {
	if(svi->point() == sp) {
	  sv = svi;
	  return true;
	}
      }
      return false;
    }

  public:
      void visit(Vertex_const_handle v) {	CGAL_NEF_TRACEN( "Vertices " << v->point() );}
      void visit(Halfedge_const_handle ) {}
      void visit(SHalfedge_const_handle ) {}
      void visit(SHalfloop_const_handle ) {}
      void visit(SFace_const_handle sf) {
	
	typename Nef_polyhedron_3::SHalfedge_const_handle sec = sf->sface_cycles_begin();
	
	int circles = 1;
	Sphere_circle first, current;
	first = current = sec->circle();
	CGAL_NEF_TRACEN( "first+current:" << first << "+" << current );
	typename Nef_polyhedron_3::SHalfedge_around_sface_const_circulator sfc(sec), send(sfc);
	CGAL_For_all(sfc, send) {
	  CGAL_NEF_TRACEN( "sedge->cirlce() " << sfc->circle() );
	  if(sfc->circle() != current) {
	    if(sfc->circle() != first)
	      ++circles;
	    current = sfc->circle();
	  }
	}
	
	CGAL_NEF_TRACEN( "first+current:" << first << "+" << current );
	CGAL_NEF_TRACEN( "circles " << circles );

	if(circles < 3)
	  omit_vertex[sf->center_vertex()] = true;
      }
  
      void visit(Halffacet_const_handle f) {

	CGAL_NEF_TRACEN( "SVertex_creator " << f->twin()->plane() );

	Halffacet_cycle_const_iterator fc = f->twin()->facet_cycles_begin();
	SHalfedge_const_handle se(fc);
	SHalfedge_around_facet_const_circulator hc(se), hend(hc);

	bool verge_found(false);
	CGAL_For_all(hc,hend) {
	  Point_3 p(hc->source()->source()->point());
	  CGAL_NEF_TRACEN(" hc " << CGAL::to_double(p.hx()) << 
			  " " << CGAL::to_double(p.hy()) <<
			  " " << CGAL::to_double(p.hz()));
	  CGAL_NEF_TRACEN(" hc->snext()->circle() " << normalized(hc->snext()->circle()));
	  if(normalized(hc->snext()->circle()) == normalized(hc->circle())) {
	    next[hc] = hc->snext();
	    CGAL_NEF_TRACEN( "set next " << hc->source()->source()->point() << ":"
		      << hc->source()->point() << "->" << hc->twin()->source()->point() << " | " 
		      << hc->snext()->source()->point() << "->" << hc->snext()->twin()->source()->point() );
 	    omit_edge[hc->twin()->source()] = omit_edge[hc->twin()->source()->twin()] = true;
	  } else
	    verge_found = true;
	}

	if(!verge_found) {
	  omit_facet[f] = true;
	  return;
	}

	SVertex_handle sv;
	if(!svertex_exists(f->twin()->plane().orthogonal_vector(), sv)) {
	  sv = SM.new_svertex(f->twin()->plane().orthogonal_vector());
	  sv->point() = normalized(sv->point());
	  sv->mark() = Mark(Point_3(0,0,0), f->mark());
	} else {
	  omit_facet[f] = true;
	  if(sv->mark().boolean() && !f->mark())
	    sv->mark().set_boolean(false);
          CGAL_NEF_TRACEN("omit facet " << f->plane());
	}
	Facet2SVertex[f] = sv;
      }
  };

  template<typename Nef_polyhedron_3>
    class SEdge_creator2 {

    typedef typename Nef_polyhedron_3::Vertex_const_handle      Vertex_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle   Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle    Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle   SHalfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfloop_const_handle   SHalfloop_const_handle;
    typedef typename Nef_polyhedron_3::SFace_const_handle       SFace_const_handle;

    typedef typename SM_decorator::SHalfedge_around_svertex_circulator
      SHalfedge_around_svertex_circulator;

    typedef typename Nef_polyhedron_3::Halffacet_cycle_const_iterator   
      Halffacet_cycle_const_iterator;
    typedef typename Nef_polyhedron_3::SHalfedge_around_facet_const_circulator
      SHalfedge_around_facet_const_circulator;

    typedef CGAL::Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge_hash;
    typedef CGAL::Unique_hash_map<Halffacet_const_handle, SVertex_handle> Facet2SVertex_hash;

    SM_decorator SM;
    Edge2SEdge_hash& Edge2SEdge;
    Facet2SVertex_hash& Facet2SVertex;

  public:
    SEdge_creator2(Sphere_map* smap, Edge2SEdge_hash& E2SE, Facet2SVertex_hash& F2SV)
      : SM(smap), Edge2SEdge(E2SE), Facet2SVertex(F2SV) {}
      
  public:
    void visit(Vertex_const_handle v) {}
    void visit(Halfedge_const_handle e) {}
    void visit(SHalfedge_const_handle se) {}
    void visit(SHalfloop_const_handle sl) {}
    void visit(SFace_const_handle sf) {}
    
    void visit(Halffacet_const_handle f) {
      
      CGAL_NEF_TRACEN( "SEdge_creator2 " << f->twin()->plane() );
      Halffacet_cycle_const_iterator fc = f->twin()->facet_cycles_begin();
      SHalfedge_const_handle sef(fc);
      SHalfedge_around_facet_const_circulator hc(sef), hend(hc);
      
      CGAL_For_all(hc, hend) {
	if(hc->sprev()->facet()->plane() == 
	   hc->facet()->plane())
	  continue;
	SHalfedge_handle thetwin;
	Halfedge_const_handle e(hc->source()->twin());
	SHalfedge_handle set = Edge2SEdge[e];
	if(set == SHalfedge_handle()) {
	  thetwin = SM.new_shalfedge_pair_at_source(Facet2SVertex[f],1);
	  CGAL_NEF_TRACEN("add stub " << Facet2SVertex[f]->point()
			  << "->" << hc->twin()->snext()->facet()->plane().orthogonal_vector());
	} else {	
	  SM.link_as_target_and_append(Facet2SVertex[f], set);
	  set->mark() = set->twin()->mark() = Mark(Point_3(0,0,0), e->mark());
	  set->circle() = Sphere_circle(set->source()->point(), set->twin()->source()->point());
	  set->twin()->circle() = set->circle().opposite();
	  thetwin = set->twin();
	  CGAL_NEF_TRACEN("complete edge " << set->source()->point()
			  << "->" << set->twin()->source()->point());
	}
	Edge2SEdge[hc->source()] = thetwin;
      }
    }
  };

  template<typename Nef_polyhedron_3>
    class SEdge_creator {

    typedef typename Nef_polyhedron_3::Vertex_const_handle      Vertex_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle   Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle    Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle   SHalfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfloop_const_handle   SHalfloop_const_handle;
    typedef typename Nef_polyhedron_3::SFace_const_handle       SFace_const_handle;

    typedef typename SM_decorator::SHalfedge_around_svertex_circulator
      SHalfedge_around_svertex_circulator;

    typedef typename Nef_polyhedron_3::Halffacet_cycle_const_iterator   
      Halffacet_cycle_const_iterator;
    typedef typename Nef_polyhedron_3::SHalfedge_around_facet_const_circulator
      SHalfedge_around_facet_const_circulator;

    typedef typename Nef_polyhedron_3::Sphere_point Sphere_point;

    typedef CGAL::Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge_hash;
    typedef CGAL::Unique_hash_map<Halffacet_const_handle, SVertex_handle> Facet2SVertex_hash;
    typedef CGAL::Unique_hash_map<Halffacet_const_handle, bool> Facet2bool_hash;
    typedef CGAL::Unique_hash_map<Halfedge_const_handle, bool> Edge2bool_hash;
    typedef CGAL::Unique_hash_map<SHalfedge_const_handle, SHalfedge_const_handle> SEdge2SEdge_hash;

    SM_decorator SM;
    Edge2SEdge_hash& Edge2SEdge;
    Facet2SVertex_hash& Facet2SVertex;
    SEdge2SEdge_hash& next;
    Facet2bool_hash& omit_facet;
    Edge2bool_hash& omit_edge;

  public:
    SEdge_creator(Sphere_map* smap, Edge2SEdge_hash& E2SE, Facet2SVertex_hash& F2SV, SEdge2SEdge_hash& SE2SE, Facet2bool_hash&  F2b, Edge2bool_hash& E2b)
      : SM(smap), Edge2SEdge(E2SE), Facet2SVertex(F2SV), next(SE2SE), omit_facet(F2b), omit_edge(E2b) {}
      
  public:
    void visit(Vertex_const_handle ) {}
    void visit(Halfedge_const_handle ) {}
    void visit(SHalfedge_const_handle ) {}
    void visit(SHalfloop_const_handle ) {}
    void visit(SFace_const_handle ) {}
    
    void visit(Halffacet_const_handle f) {

      if(omit_facet[f]) {
	CGAL_NEF_TRACEN( "omit facet " << 
			 normalized(Sphere_point(f->twin()->plane().orthogonal_vector())));
	return;
      }

      CGAL_NEF_TRACEN( "SEdge_creator " << 
		       normalized(Sphere_point(f->twin()->plane().orthogonal_vector())));
      Halffacet_cycle_const_iterator fc = f->twin()->facet_cycles_begin();
      SHalfedge_const_handle se(fc);
      SHalfedge_around_facet_const_circulator hc(se), hend(hc);
      
      
      while(hc->snext()->circle() == hc->circle()) ++hc;
      CGAL_NEF_TRACEN("verge " << hc->circle() << 
		      " " << hc->snext()->circle());
      // Found first edge such that next is not defined
      Halfedge_const_handle last = hc->twin()->source();      
      do {
	++hc;
	while(next[hc]!=SHalfedge_const_handle()) hc = next[hc];
      } while(normalized(hc->twin()->source()->point()) == 
	      normalized(last->point()));
      last = hc->twin()->source();
      // now last has a good value
      Halfedge_const_handle etwin;
      do {
	++hc;
	etwin = hc->source();
	while(next[hc]!=SHalfedge_const_handle()) hc = next[hc];
      } while(normalized(hc->twin()->source()->point()) == 
	      normalized(last->point()));      
      // hc is now one interesting corner further than last
      hend = hc;
      

      SHalfedge_handle thetwin;
      do {
	Halfedge_const_handle e = hc->twin()->source();
	CGAL_NEF_TRACEN(" check next " << &*next[hc]);
	CGAL_NEF_TRACEN(" check plane " << hc->facet()->plane());
	if(normalized(e->point()) != normalized(last->point())) {
	  Edge2SEdge[etwin] = thetwin;
	  SHalfedge_handle set = Edge2SEdge[e];
	  if(set == SHalfedge_handle()) {
	    thetwin = SM.new_shalfedge_pair_at_source(Facet2SVertex[f],1);
            CGAL_NEF_TRACEN("add stub " << Facet2SVertex[f]->point()
                            << "->" << 
			    normalized(Sphere_point(hc->snext()->facet()->plane().orthogonal_vector())));
	  } else {	
	    SM.link_as_target_and_append(Facet2SVertex[f], set);
	    set->mark() = set->twin()->mark() = Mark(Point_3(0,0,0), e->mark());
	    set->circle() = Sphere_circle(set->source()->point(), set->twin()->source()->point());
	    set->twin()->circle() = set->circle().opposite();
	    thetwin = set->twin();
            CGAL_NEF_TRACEN("complete edge " << set->source()->point()
                            << "->" << set->twin()->source()->point());
	  } 
	  last = hc->twin()->source();
	} else CGAL_NEF_TRACEN( "omit " );

	++hc;
	etwin = hc->source();
	while(next[hc]!=SHalfedge_const_handle()) hc = next[hc];
      } while(hc!=hend);
      Edge2SEdge[etwin] = thetwin;


      /*
      do {
	if(!omit_edge[hc->source()]) {
	  CGAL_NEF_TRACEN( "edge " << hc->source()->source()->point() 
		    << ":" << hc->source()->point() );

	  Halfedge_const_handle e = hc->source();
	  SHalfedge_handle set = Edge2SEdge[e->twin()];
	  if(set == SHalfedge_handle())
	    Edge2SEdge[e] = SM.new_shalfedge_pair_at_source(Facet2SVertex[f],1);
	  else {
	    SM.link_as_target_and_append(Facet2SVertex[f], set);
	    set->mark() = set->twin()->mark() = Mark(Point_3(0,0,0));
	    set->circle() = Sphere_circle(set->source()->point(), set->twin()->source()->point());
	    set->twin()->circle() = set->circle().opposite();
	    Edge2SEdge[e] = set->twin();
	  }
	} else {
	  CGAL_NEF_TRACEN( "omit edge " << hc->source()->source()->point() 
		    << ":" << hc->source()->point() );
	}
	++hc;
      } while(hc != hend);
      */

    }
  };
  
  template<typename Nef_polyhedron_3>
    class SFace_creator2 {     
  
    typedef typename Nef_polyhedron_3::Vertex_const_handle      Vertex_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle   Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle    Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle   SHalfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfloop_const_handle   SHalfloop_const_handle;
    typedef typename Nef_polyhedron_3::SFace_const_handle       SFace_const_handle;
    
    typedef typename Nef_polyhedron_3::SFace_cycle_const_iterator SFace_cycle_const_iterator;

    typedef CGAL::Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge_hash;

    SM_decorator SM;
    Edge2SEdge_hash& Edge2SEdge;

  public:
    SFace_creator2(Sphere_map* smap, Edge2SEdge_hash& E2SE) : 
      SM(smap), Edge2SEdge(E2SE) {}

      void visit(Halfedge_const_handle ) {}
      void visit(SHalfedge_const_handle ) {}
      void visit(SHalfloop_const_handle ) {}
      void visit(Halffacet_const_handle ) {}
      void visit(Vertex_const_handle ) {}

      void visit(SFace_const_handle sf) {
	CGAL_NEF_TRACEN( "SFace_creator2 " << sf->center_vertex()->point() );

	SFace_cycle_const_iterator sfc(sf->sface_cycles_begin());
	CGAL_assertion(sfc.is_shalfedge());
   	SHalfedge_const_handle sef(sfc);
	Halfedge_const_handle e(sef->source());
	SHalfedge_around_sface_circulator 
	  sec(Edge2SEdge[e]), send(sec);
	CGAL_For_all(sec, send)
	  if(sec->source()->point() == 
	     sec->twin()->source()->point())
	    return;
	SFace_handle sf_new = SM.new_sface();
	sf_new->mark() = Mark(sf->center_vertex()->point(),
			      sf->center_vertex()->mark());
	SM.link_as_face_cycle(sec,sf_new);
      }
  };

  template<typename Nef_polyhedron_3>
    class SFace_creator {     
  
    typedef typename Nef_polyhedron_3::Vertex_const_handle      Vertex_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle   Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle    Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle   SHalfedge_const_handle;
    typedef typename Nef_polyhedron_3::SHalfloop_const_handle   SHalfloop_const_handle;
    typedef typename Nef_polyhedron_3::SFace_const_handle       SFace_const_handle;

    typedef CGAL::Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge_hash;
    typedef CGAL::Unique_hash_map<Vertex_const_handle, bool> Vertex2bool_hash;
    typedef CGAL::Unique_hash_map<SFace_const_handle, bool> SFace2bool_hash;

    const Nef_polyhedron_3& N3;
    SM_decorator SM;
    Edge2SEdge_hash& Edge2SEdge;
    Vertex2bool_hash& omit_vertex;
    SFace2bool_hash& Shell;

  public:
    SFace_creator(const Nef_polyhedron_3& N, Sphere_map* smap, Edge2SEdge_hash& E2SE, Vertex2bool_hash& V2b, SFace2bool_hash SHELL) : 
      N3(N), SM(smap), Edge2SEdge(E2SE), omit_vertex(V2b), Shell(SHELL) {}

      void visit(Halfedge_const_handle ) {}
      void visit(SHalfedge_const_handle ) {}
      void visit(SHalfloop_const_handle ) {}
      void visit(SFace_const_handle ) {}
      void visit(Halffacet_const_handle ) {}

      void visit(Vertex_const_handle v) {

	CGAL_NEF_TRACEN( "SFace_creator " << v->point() );

	if(omit_vertex[v]) {
	  CGAL_NEF_TRACEN("omit " << v->point() );
	  return;
	}

        typename Nef_polyhedron_3::Nef_polyhedron_S2 SD(N3.get_sphere_map(v));
   
        /*
	typename Nef_polyhedron_3::SFace_const_iterator sf = SD.sfaces_begin();
	while(sf != SD.sfaces_end() && !Shell[sf]) ++sf;
        CGAL_assertion(sf != SD.sfaces_end());
        */

	typename Nef_polyhedron_3::Halfedge_const_iterator ei(SD.svertices_begin());
	SHalfedge_handle se = Edge2SEdge[ei];
	while(se == SHalfedge_handle()) {
	  ++ei;
	  se = Edge2SEdge[ei];
	}

	CGAL_assertion(ei != SD.svertices_end());

	SFace_handle sf_new = SM.new_sface();
	sf_new->mark() = Mark(v->point(), v->mark());
	SM.link_as_face_cycle(se,sf_new);
      }
  };

  struct VECTOR_ADDITION {
    Mark operator()(const Mark& b1, const Mark& b2) const {
      return b1+b2;
    }
  };

  Object_handle top;
  Object_handle bottom;

  void locate_top_and_bottom() {
    std::vector<SFace_iterator> topSF;
    std::vector<SFace_iterator> bottomSF;
    SFace_iterator sfi = this->sfaces_begin();
    topSF.push_back(sfi);
    bottomSF.push_back(sfi);

    Comparison_result cr;
    for(++sfi;sfi != this->sfaces_end(); ++sfi) {
      cr = compare_z(sfi->mark(), (*topSF.begin())->mark());
      if(cr != CGAL::SMALLER) {
	if(cr == CGAL::LARGER)
	  topSF.clear();
	topSF.push_back(sfi);	
      }
      cr = compare_z(sfi->mark(), (*bottomSF.begin())->mark());
      if(cr != CGAL::LARGER) {
	if(cr == CGAL::SMALLER)
	  bottomSF.clear();
	bottomSF.push_back(sfi);	
      }    
    }

    SFace_handle sf(topSF.front());
    if(topSF.size()==1)
      top = Object_handle(SFace_const_handle(sf));
    else {
      SHalfedge_handle se(sf->sface_cycles_begin());
      SHalfedge_around_sface_circulator sfc(se), sfend(sfc);       
      
      if(topSF.size()==2) {
	while(sfc->circle().c()!=0) {
	  ++sfc;
	  CGAL_assertion(sfc != sfend);
	}
	top = Object_handle(SHalfedge_const_handle(sfc));
      } else {
	CGAL_assertion(topSF.size() > 0);
	while(sfc->source()->point().hx()!=0 || 
	      sfc->source()->point().hy()!=0 ||
	      sfc->source()->point().hz()<=0) {
	  ++sfc;
	  CGAL_assertion(sfc != sfend);
	}
	top = Object_handle(SVertex_const_handle(sfc->source()));      
      }
    }

    sf = bottomSF.front();
    if(bottomSF.size()==1)
      bottom = Object_handle(SFace_const_handle(sf));
    else {
      SHalfedge_handle se(sf->sface_cycles_begin());
      SHalfedge_around_sface_circulator sfc(se), sfend(sfc);       
      
      if(bottomSF.size()==2) {
	while(sfc->circle().c()!=0) {
	  ++sfc;
	  CGAL_assertion(sfc != sfend);
	}
	bottom = Object_handle(SHalfedge_const_handle(sfc));
      } else {
	CGAL_assertion(bottomSF.size() > 0);
	while(sfc->source()->point().hx()!=0 || 
	      sfc->source()->point().hy()!=0 || 
	      sfc->source()->point().hz()>=0) {
	  ++sfc;
	  CGAL_assertion(sfc != sfend);
	}
	bottom = Object_handle(SVertex_const_handle(sfc->source()));      
      }
    }
  }

  void erase_redundant_vertices() {

    std::cerr << "erase redundant vertices " << std::endl;

    /*
    SVertex_iterator svi;
    CGAL_forall_svertices(svi, *this) {
      bool erase(true);
      SHalfedge_around_svertex_circulator 
	svc(svi->out_sedge()), send(svc);
      CGAL_For_all(svc, send) {
	if(svc->incident_sface() != SFace_handle() ||
	   svc->twin()->incident_sface() != SFace_handle()) {
	  erase = false;
	  break;
	}
      }
      if(!erase) continue;
      SHalfedge_handle se(svi->out_sedge());
      while(se->twin()->snext() != se) {
	SHalfedge_handle se_next(se->twin()->snext());
	delete_edge_pair(se);
	se = se_next;
      }
      delete_edge_pair(se);
      delete_vertex_only(svi);
    }
    */

    std::list<SHalfedge_handle> redundant;
    SHalfedge_iterator sei;
    CGAL_forall_sedges(sei, *this) {
      if(sei->source()->point() ==
	 sei->twin()->source()->point()) {
	redundant.push_back(sei);
      }
    }

    CGAL::Unique_hash_map<SHalfedge_handle, bool> erased(false);
    typename std::list<SHalfedge_handle>::iterator ri;
    for(ri = redundant.begin(); ri != redundant.end(); ++ri) {
      if(erased[(*ri)]) continue;

      SVertex_handle src((*ri)->source());
      SVertex_handle tgt((*ri)->twin()->source());

      std::cerr << "erase " << src->point() << std::endl;
      std::cerr << &*src << " " << &*tgt << std::endl;

      SHalfedge_handle prev((*ri)->sprev());
      SHalfedge_handle next((*ri)->snext());

      std::cerr << "prev " << &*prev->source()
		<< "->" << &*prev->twin()->source() << std::endl;      
      std::cerr << "next " << &*next->source()
		<< "->" << &*next->twin()->source() << std::endl;
 
      if(prev->source() ==
	 next->twin()->source()) {
	std::cerr << "delete before" << std::endl;
	SHalfedge_handle sein;
	SHalfedge_handle se_cas(next->twin()->snext());
	SFace_handle sf(next->twin()->incident_sface());
	if(sf != SFace_handle()) {
	  std::cerr << "not null " << std::endl;
	  SFace_cycle_iterator sfci(sf->sface_cycles_begin());
	  CGAL_assertion(sfci.is_shalfedge());
	  sein = sfci;
	  undo_sm_boundary_object(sein, sf);
	  sein = next->twin()->sprev();
	  std::cerr << "sein " << sein->source()->point()
		    << "->" << sein->twin()->source()->point()
		    << std::endl;
	}
	erased[next] = erased[next->twin()] = true;
	delete_edge_pair(next);
	next = se_cas;
	if(sf != SFace_handle()) {
	  link_as_face_cycle(sein, sf);
	}
      }
 
      prev->snext() = next;
      next->sprev() = prev;

      SHalfedge_handle tprev((*ri)->twin()->sprev());
      SHalfedge_handle tnext((*ri)->twin()->snext());
      if(tprev->source() ==
	 tnext->twin()->source()) {
	std::cerr << "delete twin" << std::endl;
	SHalfedge_handle sein;
	SHalfedge_handle se_cap(tprev->twin()->sprev());
	SFace_handle sf(tprev->twin()->incident_sface());
	if(sf != SFace_handle()) {
	  SFace_cycle_iterator sfci(sf->sface_cycles_begin());
	  CGAL_assertion(sfci.is_shalfedge());
	  sein=sfci;
	  undo_sm_boundary_object(sein, sf);
	  sein = tprev->twin()->snext();
	}
	erased[tprev] = erased[tprev->twin()] = true;
	delete_edge_pair(tprev);
	tprev = se_cap;
	if(sf != SFace_handle())
	  link_as_face_cycle(sein, sf);
      }
      tprev->snext() = tnext;
      tnext->sprev() = tprev;

      std::cerr << "tprev " << &*tprev->source()
		<< "->" << &*tprev->twin()->source() << std::endl;      
      std::cerr << "tnext " << &*tnext->source()
		<< "->" << &*tnext->twin()->source() << std::endl;

      while(next != tnext) {
      std::cerr << "next " << &*next->source()
		<< "->" << &*next->twin()->source() << std::endl;
	SHalfedge_handle se_cas(next->twin()->snext());
	next->source() = src;
	next = se_cas;
      }
      delete_edge_pair_only(*ri);
      delete_vertex_only(tgt);
    }
  }

 public:
  Gaussian_map() : Base(new Sphere_map) {}

  template<typename NK, typename Items> 
  Gaussian_map(const CGAL::Nef_polyhedron_3<NK, Items>& N3,
	      typename CGAL::Nef_polyhedron_3<NK, Items>::Volume_const_iterator c) : Base(new Sphere_map) {

    typedef CGAL::Nef_polyhedron_3<NK, Items> Nef_polyhedron_3;
    typedef typename Nef_polyhedron_3::Vertex_const_iterator 
      Vertex_const_iterator;
    typedef typename Nef_polyhedron_3::Halffacet_const_iterator
      Halffacet_const_iterator;
    typedef typename Nef_polyhedron_3::Halffacet_cycle_const_iterator
      Halffacet_cycle_const_iterator;
    typedef typename Nef_polyhedron_3::SHalfedge_around_facet_const_circulator
      SHalfedge_around_facet_const_circulator;
    typedef typename Nef_polyhedron_3::Vertex_const_handle
      Vertex_const_handle;
    typedef typename Nef_polyhedron_3::Volume_const_handle
      Volume_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle
      Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle
      Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle
      SHalfedge_const_handle;
    typedef typename Nef_polyhedron_3::SFace_const_handle
      SFace_const_handle;

    Unique_hash_map<Halffacet_const_handle, SVertex_handle> Facet2SVertex;
    Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge;
    Unique_hash_map<Halffacet_const_handle, bool> Facet2bool;
    Unique_hash_map<Vertex_const_handle, bool> Vertex2bool(false);
    Unique_hash_map<Halfedge_const_handle, bool> Edge2bool(false);
    Unique_hash_map<SHalfedge_const_handle, SHalfedge_const_handle> SEdge2SEdge;
    Unique_hash_map<SFace_const_handle, bool> Shell(false);

    SFace_const_handle sf = c->shells_begin();

    SVertex_creator<Nef_polyhedron_3> create_svertices(this->sphere_map(), Facet2SVertex, Facet2bool, SEdge2SEdge, Vertex2bool, Edge2bool, Shell);
    SEdge_creator<Nef_polyhedron_3>   create_sedges(this->sphere_map(), Edge2SEdge, Facet2SVertex, SEdge2SEdge, Facet2bool, Edge2bool);
    SFace_creator<Nef_polyhedron_3>   create_sfaces(N3, this->sphere_map(), Edge2SEdge, Vertex2bool, Shell);

    N3.visit_shell_objects(sf, create_svertices);
    N3.visit_shell_objects(sf, create_sedges);
    N3.visit_shell_objects(sf, create_sfaces);
  }

  template<typename NK> 
    Gaussian_map(const CGAL::Nef_polyhedron_3<NK>& N3) : Base(new Sphere_map) {

    typedef CGAL::Nef_polyhedron_3<NK> Nef_polyhedron_3;
    typedef typename Nef_polyhedron_3::Vertex_const_iterator 
      Vertex_const_iterator;
    typedef typename Nef_polyhedron_3::Halffacet_const_iterator
      Halffacet_const_iterator;
    typedef typename Nef_polyhedron_3::Halffacet_cycle_const_iterator
      Halffacet_cycle_const_iterator;
    typedef typename Nef_polyhedron_3::SHalfedge_around_facet_const_circulator
      SHalfedge_around_facet_const_circulator;
    typedef typename Nef_polyhedron_3::Volume_const_handle
      Volume_const_handle;
    typedef typename Nef_polyhedron_3::Halfedge_const_handle
      Halfedge_const_handle;
    typedef typename Nef_polyhedron_3::Halffacet_const_handle
      Halffacet_const_handle;
    typedef typename Nef_polyhedron_3::SHalfedge_const_handle
      SHalfedge_const_handle;
    
    Unique_hash_map<Halffacet_const_handle, SVertex_handle> Facet2SVertex;
    Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge;

    Volume_const_handle c(--N3.volumes_end());

    Halffacet_const_iterator f;
    CGAL_forall_halffacets(f,N3) {
      if(f->incident_volume() != c) continue;
      SVertex_handle sv = new_svertex(f->twin()->plane().orthogonal_vector());
      sv->mark() = Mark(Point_3(0,0,0), f->mark());
      Facet2SVertex[f] = sv;
    }
    
    CGAL_forall_halffacets(f,N3) {
      if(f->incident_volume() != c) continue;
      Halffacet_cycle_const_iterator fc = f->twin()->facet_cycles_begin();
      SHalfedge_const_handle se(fc);
      SHalfedge_around_facet_const_circulator hc(se), hend(hc);
      do {
	Halfedge_const_handle e = hc->source();
	SHalfedge_handle set = Edge2SEdge[e->twin()];
	if(set == SHalfedge_handle())
	  Edge2SEdge[e] = new_shalfedge_pair_at_source(Facet2SVertex[f],1);
	else {
	  link_as_target_and_append(Facet2SVertex[f], set);
	  set->mark() = set->twin()->mark() = Mark(Point_3(0,0,0), e->mark());
	  set->circle() = Sphere_circle(set->source()->point(), set->twin()->source()->point());
	  set->twin()->circle() = set->circle().opposite();
	  Edge2SEdge[e] = set->twin();
	}
	++hc;
      } while(hc != hend);
    }
    
    Vertex_const_iterator v;
    CGAL_forall_vertices(v,N3) {

      typename Nef_polyhedron_3::Nef_polyhedron_S2 SD(N3.get_sphere_map(v));
      Halfedge_const_handle e(SD.svertices_begin());
      SHalfedge_handle se = Edge2SEdge[e];
      SFace_handle sf = this->new_sface();
      sf->mark() = Mark(v->point(), v->mark());
      link_as_face_cycle(se,sf);
    }
  }
  
  template<typename PK> 
  Gaussian_map(const CGAL::Polyhedron_3<PK>& P, bool closed = true) 
    : Base(new Sphere_map) {
    
    typedef CGAL::Polyhedron_3<PK> Polyhedron_3;
    typedef typename Polyhedron_3::Vertex_const_iterator 
      Vertex_const_iterator;
    typedef typename Polyhedron_3::Facet_const_iterator
      Facet_const_iterator;
    typedef typename Polyhedron_3::Halfedge_around_facet_const_circulator
      Halfedge_around_facet_const_circulator;
    typedef typename Polyhedron_3::Halfedge_const_handle
      Halfedge_const_handle;
    typedef typename Polyhedron_3::Facet_const_handle
      Facet_const_handle;
    
    Unique_hash_map<Facet_const_handle, SVertex_handle> Facet2SVertex;
    Unique_hash_map<Halfedge_const_handle, SHalfedge_handle> Edge2SEdge;

    Facet_const_iterator f;
    for(f = P.facets_begin(); f != P.facets_end(); ++f) {
      SVertex_handle sv = new_svertex(f->plane().orthogonal_vector());
      sv->mark() = Mark(Point_3(0,0,0), closed);
      Facet2SVertex[f] = sv;
    }

    for(f = P.facets_begin(); f != P.facets_end(); ++f) {
      Halfedge_around_facet_const_circulator hc(f->facet_begin()),hend(hc);
      do {
	Halfedge_const_handle e = hc;
	SHalfedge_handle set = Edge2SEdge[e->opposite()];
	if(set == SHalfedge_handle())
	  Edge2SEdge[e] = new_shalfedge_pair_at_source(Facet2SVertex[f],1);
	else {
	  link_as_target_and_append(Facet2SVertex[f], set,1);
	  set->mark() = set->twin()->mark() = Mark(Point_3(0,0,0), closed);
	  set->circle() = Sphere_circle(set->source()->point(), set->twin()->source()->point());
	  set->twin()->circle() = set->circle().opposite();
	  Edge2SEdge[e] = set->twin();
	}
	++hc;
      } while(hc != hend);
    }
   
    Vertex_const_iterator v;
    for(v = P.vertices_begin(); v != P.vertices_end(); ++v) {
      Halfedge_const_handle e(v->halfedge());
      SHalfedge_handle se = Edge2SEdge[e];
      SFace_handle sf = this->new_sface();
      sf->mark() = Mark(v->point(), closed);
      link_as_face_cycle(se,sf);
    }
  }

  Gaussian_map(typename Nef::Halffacet_const_handle f) 
    : Base(new Sphere_map) 
  {
    SVertex_handle 
      sv1(this->new_svertex(CGAL::ORIGIN + 
			    f->plane().orthogonal_vector())),
      sv2(this->new_svertex(CGAL::ORIGIN-sv1->point()));
	  
    sv1->mark() = sv1->mark() = 
      Mark(Point_3(0,0,0), f->mark());

    typename Nef::SHalfedge_around_facet_const_circulator 
      sfc(f->facet_cycles_begin()), send(sfc);

    SHalfedge_handle se1 = this->new_shalfedge_pair(sv1, sv2);
    Point_3 orth(sfc->twin()->source()->point());
    se1->circle() = Sphere_circle(orth.hx(), 
				  orth.hy(),
				  orth.hz());
    se1->twin()->circle() = se1->circle().opposite();
    se1->mark() = se1->twin()->mark() = 
      Mark(Point_3(0,0,0), sfc->twin()->source()->mark());
    
    SHalfedge_handle se_prev(se1);

    ++sfc;
    CGAL_For_all(sfc, send) {
      SHalfedge_handle 
	se(this->new_shalfedge_pair(se_prev, se_prev->twin(),
				    SM_decorator::AFTER, SM_decorator::BEFORE));
      Point_3 orth(sfc->twin()->source()->point());
      se->circle() = Sphere_circle(orth.hx(), 
				   orth.hy(),
				   orth.hz());
      se->twin()->circle() = se->circle().opposite();
      se->mark() = se->twin()->mark() = 
	Mark(Point_3(0,0,0), sfc->twin()->source()->mark());      
      se_prev = se;
    }

    typename Sphere_map::SHalfedge_around_svertex_circulator 
      svc(sv1->out_sedge());
    CGAL_For_all(sfc, send) {
      SFace_handle sf = this->new_sface();
      sf->mark() = 
	Mark(sfc->source()->source()->point(), sfc->source()->source()->mark());
      this->link_as_face_cycle(svc->twin(), sf);
      ++svc;
    } 
  }

  Gaussian_map(typename Nef::Halfedge_const_handle e) 
    : Base(new Sphere_map) 
  {
    SHalfloop_handle sl = this->new_shalfloop_pair();
    Point_3 p = e->twin()->point();
    sl->circle() = Sphere_circle(p.hx(), p.hy(), p.hz());
    sl->twin()->circle() = sl->circle().opposite();
    sl->mark() = sl->twin()->mark() = 
      Mark(Point_3(0,0,0), e->mark());
    
    SFace_handle sf1 = this->new_sface();
    SFace_handle sf2 = this->new_sface();
    sf1->mark() = Mark(e->source()->point(), 
		       e->source()->mark());
    sf2->mark() = Mark(e->twin()->source()->point(),
		       e->twin()->source()->mark());

    this->link_as_loop(sl, sf1);
    this->link_as_loop(sl->twin(), sf2);
  }

  Gaussian_map(typename Nef::Vertex_const_handle v) 
    : Base(new Sphere_map) 
  {
    SFace_handle sf = this->new_sface();
    sf->mark() = Mark(v->point(), v->mark());
  }


  void simplify() 
  {
	  CGAL_NEF_TRACEN("simplify");
    
	  typedef typename CGAL::Union_find<SFace_handle>::handle Union_find_handle;
	  CGAL::Unique_hash_map< SFace_handle, Union_find_handle> Pitem(NULL);
	  CGAL::Union_find< SFace_handle> UF;
    
	  SFace_iterator f;
	  CGAL_forall_sfaces(f,*this) {
		  Pitem[f] = UF.make_set(f);
		  clear_face_cycle_entries(f);
	  }
    
	  SHalfedge_iterator e, en;
	  for(e = this->shalfedges_begin(); e != this->shalfedges_end(); e = en) 
	  { 
		  en = e; ++en; if ( en==e->twin() ) ++en;
		  CGAL_NEF_TRACEN("can simplify ? " << PH(e));
		  CGAL_NEF_TRACEN(e->mark() << " " << e->incident_sface()->mark() 
						  << " " << e->twin()->incident_sface()->mark());
		  if (e->incident_sface()->mark() == 
			  e->twin()->incident_sface()->mark()) {
			  CGAL_NEF_TRACEN("deleting "<<PH(e));
			  if ( !UF.same_set(Pitem[e->incident_sface()],
								Pitem[e->twin()->incident_sface()]) ) {
	  
				  UF.unify_sets( Pitem[e->incident_sface()],
								 Pitem[e->twin()->incident_sface()] );
				  CGAL_NEF_TRACEN("unioning disjoint faces");
			  }
	
			  CGAL_NEF_TRACEN("is_closed_at_source " << is_closed_at_source(e) << 
							  " " << is_closed_at_source(e->twin()));
			  delete_edge_pair(e);
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
	  for(v = this->svertices_begin(); v != this->svertices_end(); v=vn) 
	  {
		  vn=v; ++vn;
		  if ( is_isolated(v) ) {
			  delete_vertex_only(v);
			  continue;
		  }
		  if ( has_outdeg_two(v)) {
			  merge_edge_pairs_at_target(first_out_edge(v)->sprev()); 
		  } 
	  }
    
	  SFace_iterator fn;
	  for (f = fn = this->sfaces_begin(); f != this->sfaces_end(); f=fn) 
	  { 
		  ++fn;
		  Union_find_handle pit = Pitem[f];
		  if ( UF.find(pit) != pit ) {
			  CGAL_NEF_TRACEN("delete face " << &*f);
			  delete_face_only(f);
		  }
	  }
  }      
  
  void minkowski_sum(const Gaussian_map& G1, const Gaussian_map& G2) {
    SM_overlayer O(this->sphere_map());
#ifdef CGAL_NEF3_TIMER_OVERLAY
    CGAL::Timer t;
    t.start();
#endif // CGAL_NEF3_TIMER_OVERLAY
    O.subdivide(G1.sphere_map(), G2.sphere_map(), true);
#ifdef CGAL_NEF3_TIMER_OVERLAY
    t.stop();
    std::cout << "Runtime_overlay " << t.time() << std::endl;
#endif // CGAL_NEF3_TIMER_OVERLAY
    VECTOR_ADDITION va;
    O.select(va);
    simplify();
  }
  
  Object_handle get_top() {
    return top;
  }
  Object_handle get_bottom() {
    return bottom;
  }
  
  Object_handle locate_top() {
    std::vector<SFace_iterator> topSF;
    SFace_iterator sfi = this->sfaces_begin();
    topSF.push_back(sfi);
    
    Comparison_result cr;
    for(++sfi;sfi != this->sfaces_end(); ++sfi) {
      cr = compare_z(sfi->mark(), (*topSF.begin())->mark());
      if(cr != CGAL::SMALLER) {
	if(cr == CGAL::LARGER)
	  topSF.clear();
	topSF.push_back(sfi);	
      }
    }
    
    SFace_handle sf(topSF.front());
    if(topSF.size()==1)
    return Object_handle(SFace_const_handle(sf));
    else {
      SHalfedge_handle se(sf->sface_cycles_begin());
      SHalfedge_around_sface_circulator sfc(se), sfend(sfc);       
      
      if(topSF.size()==2) {
	while(sfc->circle().c()!=0) {
	  ++sfc;
	  CGAL_assertion(sfc != sfend);
	}
	return Object_handle(SHalfedge_const_handle(sfc));
      } else {
	CGAL_assertion(topSF.size() > 0);
	while(sfc->source()->point().hx()!=0 || 
	      sfc->source()->point().hy()!=0 ||
	      sfc->source()->point().hz()<=0) {
	  ++sfc;
	  CGAL_assertion(sfc != sfend);
	}
	return Object_handle(SVertex_const_handle(sfc->source()));      
      }
    }
    CGAL_error_msg("line should not be executed");
    return Object_handle();
  }
  
  Object_handle locate_bottom() {
    std::vector<SFace_iterator> bottomSF;
    SFace_iterator sfi = this->sfaces_begin();
    bottomSF.push_back(sfi);
    
    Comparison_result cr;
    for(++sfi;sfi != this->sfaces_end(); ++sfi) {
      cr = compare_z(sfi->mark(), (*bottomSF.begin())->mark());
      if(cr != CGAL::LARGER) {
	if(cr == CGAL::SMALLER)
	  bottomSF.clear();
	bottomSF.push_back(sfi);	
      }
    }
    
    SFace_handle sf(bottomSF.front());
    if(bottomSF.size()==1)
    return Object_handle(SFace_const_handle(sf));
    else {
      SHalfedge_handle se(sf->sface_cycles_begin());
      SHalfedge_around_sface_circulator sfc(se), sfend(sfc);       
      
	if(bottomSF.size()==2) {
	  while(sfc->circle().c()!=0) {
	    ++sfc;
	    CGAL_assertion(sfc != sfend);
	  }
	  return Object_handle(SHalfedge_const_handle(sfc));
	} else {
	  CGAL_assertion(bottomSF.size() > 0);
	  while(sfc->source()->point().hx()!=0 || 
		sfc->source()->point().hy()!=0 ||
		sfc->source()->point().hz()>=0) {
	    ++sfc;
	    CGAL_assertion(sfc != sfend);
	  }
	  return Object_handle(SVertex_const_handle(sfc->source()));      
	}
    }
    CGAL_error_msg("line should not be executed");
    return Object_handle();
  }
  
};

template<typename Kernel, typename Nef, typename Mark>
std::ostream& operator<<(std::ostream& out, const CGAL::Gaussian_map<Kernel, Nef, Mark>& G) {
  out << "OFF" << std::endl;
  out << G.number_of_sfaces() << " " << G.number_of_svertices() << " 0" << std::endl;
  
  typedef typename CGAL::Gaussian_map<Kernel, Nef, Mark>::SFace_const_iterator 
    SFace_const_iterator;
  CGAL::Unique_hash_map<SFace_const_iterator, int> SFace2int;
  
  int i=0;
  SFace_const_iterator sf;
  CGAL_forall_sfaces(sf, G) {
    SFace2int[sf] = i++;
    out << CGAL::to_double(sf->mark().x()) << " " 
	<< CGAL::to_double(sf->mark().y()) << " " 
	<< CGAL::to_double(sf->mark().z()) << std::endl;
  }

  typename CGAL::Gaussian_map<Kernel, Nef, Mark>::SVertex_const_iterator sv;
  CGAL_forall_svertices(sv,G) {
    typename CGAL::Gaussian_map<Kernel, Nef, Mark>::SHalfedge_around_svertex_const_circulator 
      svc(G.first_out_edge(sv)),
      svc1(svc),
      svend(svc);
    out << std::distance(++svc1,svend)+1;
    CGAL_For_all(svc,svend)
      out << " " << SFace2int[svc->incident_sface()];
    out << std::endl;
  }

  return out;
}

CGAL_END_NAMESPACE
#endif // CGAL_MS3_GAUSSIAN_MAP
