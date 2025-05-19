// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_CD3_REFLEX_VERTEX_SEARCHER_H
#define CGAL_CD3_REFLEX_VERTEX_SEARCHER_H

#include <CGAL/license/Convex_decomposition_3.h>


#include<CGAL/Nef_3/SNC_decorator.h>
#include<CGAL/Convex_decomposition_3/is_reflex_sedge.h>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 229
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template<typename Nef_>
class Reflex_vertex_searcher : public Modifier_base<typename Nef_::SNC_structure> {

  typedef Nef_                                            Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_structure          SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>              SNC_decorator;

  typedef typename SNC_structure::Sphere_map Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>  SM_decorator;
  typedef CGAL::SM_point_locator<SM_decorator>    SM_point_locator;

  typedef typename SNC_structure::Vertex_handle           Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle         Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle        Halffacet_handle;
  typedef typename SNC_structure::SHalfedge_handle        SHalfedge_handle;
  typedef typename SNC_structure::SHalfloop_handle        SHalfloop_handle;
  typedef typename SNC_structure::SFace_handle            SFace_handle;

  typedef typename SNC_structure::SVertex_handle    SVertex_handle;

  typedef typename SNC_structure::Vertex_iterator         Vertex_iterator;
  typedef typename SNC_structure::Volume_iterator         Volume_iterator;
  typedef typename SNC_structure::SHalfedge_iterator      SHalfedge_iterator;

  typedef typename SNC_structure::Shell_entry_iterator    Shell_entry_iterator;
  typedef typename SNC_structure::SFace_cycle_iterator    SFace_cycle_iterator;
  typedef typename SNC_structure::SHalfedge_around_sface_circulator SHalfedge_around_sface_circulator;

  typedef typename SNC_structure::Sphere_point            Sphere_point;
  typedef typename SNC_structure::Sphere_segment          Sphere_segment;

 public:
  typedef CGAL::Unique_hash_map<Vertex_handle, int>       Reflex_vertex_map;

 private:
  struct Reflex_vertex_visitor {

    Sphere_point dir;
    Reflex_vertex_map& vertex_map;

    Reflex_vertex_visitor(Sphere_point dir_in,
                          Reflex_vertex_map& vm)
      : dir(dir_in), vertex_map(vm) {}

    void visit(Vertex_handle /*v*/) const {}
    void visit(Halfedge_handle /*e*/) const {}
    void visit(Halffacet_handle /*f*/) const {}
    void visit(SHalfloop_handle /*sl*/) const {}
    void visit(SFace_handle /*sf*/) const {}
    void visit(SHalfedge_handle se) const {
      //      if(vertex_map[se->source()->source()]==3) return;
      int isrse = is_reflex_sedge<SNC_structure>(se, dir);
      if(isrse==0) return;
      vertex_map[se->source()->source()] |= isrse;
      vertex_map[se->prev()->source()->source()] |= isrse;
    }
  };

 public:
  Sphere_point dir;
  Reflex_vertex_map vertex_map;
  Vertex_iterator pos, neg, begin, end;

  Reflex_vertex_searcher(Sphere_point dir_in)
      : dir(dir_in), vertex_map(0) {}

  /*
  void operator()(SNC_structure& snc) {
    pos = neg = begin = snc.vertices_begin();
    end = snc.vertices_end();

    Reflex_vertex_visitor rev(dir, vertex_map);
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

  int is_reflex_vertex(Vertex_handle vi) {
    int result = 0;
    SM_point_locator PL(&*vi);
    Object_handle o[2] = {PL.locate(dir), PL.locate(dir.antipode())};

    for(int i=0; i<2; ++i) {
      SFace_handle sf;
      bool markedsf = assign(sf, o[i]) && sf->mark();

      CGAL_NEF_TRACEN("markedsf " << markedsf);
      CGAL_NEF_TRACEN("sf " << &*sf);

      if(markedsf) {
        SFace_cycle_iterator sfci(sf->sface_cycles_begin());
        for(; sfci != sf->sface_cycles_end(); ++sfci) {
          SHalfedge_around_sface_circulator sfc(sfci), send(sfc);
          CGAL_For_all(sfc, send) {
            int isrse = is_reflex_sedge<SNC_structure>(sfc, dir, false);
            if(isrse == 0) continue;
            isrse &= (i + 1);
            result |= isrse;
          }
        }
      }
    }

    vertex_map[vi]=result;
    return result;
  }

  bool need_to_shoot(SVertex_handle sv, bool turn_around) {
    Sphere_point upDown = turn_around ? dir.antipode() : dir;
    Sphere_segment sray(sv->point(), upDown);
    SM_point_locator smpl(&*sv->source());
    Sphere_point ip;
    Object_handle o = smpl.ray_shoot(sray, ip, false, false, true);
    SHalfedge_handle se;
    return !(assign(sv, o) || assign(se, o));
  }

  void operator()(SNC_structure& snc) {
    pos = neg = begin = snc.vertices_begin();
    end = snc.vertices_end();

    Vertex_iterator vi;
    CGAL_forall_vertices(vi, snc)
      vertex_map[vi] |= is_reflex_vertex(vi);
  }

  void set_back_positive_reflex_vertices() {
    pos = begin;
  }

  void set_back_negative_reflex_vertices() {
    neg = begin;
  }

  bool get_next_positive_reflex_vertex(Vertex_handle& v) {
    if(pos == end) return false;
    while(pos != end && (vertex_map[pos]&1)!=1)
      ++pos;
    if(pos == end) return false;
    v = pos;
    ++pos;
    return true;
  }

  bool get_next_negative_reflex_vertex(Vertex_handle& v) {
    if(neg == end) return false;
    while(neg != end && (vertex_map[neg]&2)!=2)
      ++neg;
    if(neg == end) return false;
    v = neg;
    ++neg;
    /*
    SHalfedge_iterator sei;
    for(sei = v->shalfedges_begin();
        sei != v->shalfedges_end(); ++sei) {
      int isrse = is_reflex_sedge<SNC_structure>(sei, dir);
      std::cerr << "check " << sei->source()->source()->point()
                << "->" << sei->source()->twin()->source()->point()
                << ": " << isrse << std::endl;
    }
    */
    return true;
  }
};

} //namespace CGAL
#endif // CGAL_CD3_REFLEX_VERTEX_SEARCHER_H
