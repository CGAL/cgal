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
#ifndef CGAL_MINKOWSKI_GAUSSIAN_MAP_TO_NEF_3_H
#define CGAL_MINKOWSKI_GAUSSIAN_MAP_TO_NEF_3_H

#include <CGAL/license/Minkowski_sum_3.h>


#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Minkowski_sum_3/Gaussian_map.h>
#include <CGAL/Modifier_base.h>

namespace CGAL {

template<typename Nef3>
  class Gaussian_map_to_nef_3 : public Modifier_base<typename Nef3::SNC_structure > {

  typedef typename Nef3::Kernel                   Kernel;
  typedef typename Nef3::SNC_structure            SNC_structure;
  typedef typename SNC_structure::Sphere_map      Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>          SM_decorator;
  typedef CGAL::Gaussian_map<Kernel, Nef3>        Gaussian_map;

  typedef typename Gaussian_map::SFace_const_iterator       SFace_const_iterator;
  typedef typename Gaussian_map::SFace_const_handle         SFace_const_handle;
  typedef typename Gaussian_map::SHalfedge_const_iterator   SHalfedge_const_iterator;
  typedef typename Gaussian_map::SHalfedge_const_handle     SHalfedge_const_handle;
  typedef typename Gaussian_map::SHalfloop_const_handle     SHalfloop_const_handle;
  typedef typename Gaussian_map::SVertex_const_iterator     SVertex_const_iterator;
  typedef typename Gaussian_map::SVertex_const_handle       SVertex_const_handle;
  typedef typename Gaussian_map::SHalfedge_around_sface_const_circulator
    SHalfedge_around_sface_const_circulator;

  typedef typename SNC_structure::Vertex_handle            Vertex_handle;
  typedef typename SNC_structure::SVertex_handle           SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle         SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle             SFace_handle;

  typedef typename SNC_structure::Sphere_circle            Sphere_circle;
  typedef typename SNC_structure::Sphere_point             Sphere_point;

  const Gaussian_map& G;

 public:
  Gaussian_map_to_nef_3(const Gaussian_map& Gin) : G(Gin) {}

  void create_solid(SNC_structure& snc) {

    CGAL::Unique_hash_map<SHalfedge_const_handle, int> SE2i;
    SHalfedge_const_iterator sei;
    CGAL_forall_sedges(sei, G) {
      SE2i[sei] = Index_generator::get_unique_index();
      SE2i[sei->twin()] = SE2i[sei];
    }

    CGAL::Unique_hash_map
      <SVertex_const_handle, std::pair<int, int> > SV2i;
    SVertex_const_iterator svi;
    CGAL_forall_svertices(svi, G)
      SV2i[svi] = std::pair<int, int>
      (Index_generator::get_unique_index(),
       Index_generator::get_unique_index());

    CGAL::Unique_hash_map<SFace_const_handle, Vertex_handle> sface2vertex;
    SFace_const_iterator sfi;
    for(sfi = G.sfaces_begin(); sfi != G.sfaces_end(); ++sfi) {
      sface2vertex[sfi] = snc.new_vertex(sfi->mark().point(),
                                         sfi->mark().boolean());
    }

    for(sfi = G.sfaces_begin(); sfi != G.sfaces_end(); ++sfi) {
      Vertex_handle v = sface2vertex[sfi];
      SM_decorator SM(&*v);

      SHalfedge_const_handle sec = sfi->sface_cycles_begin();
      SHalfedge_around_sface_const_circulator sfc(sec), sfcend(sfc);

      SVertex_handle sv, sv_prev, sv_first;
      SHalfedge_handle se, se_prev, se_first;
      sv_first =
        SM.new_svertex(sface2vertex[sfc->twin()->incident_sface()]->point()-v->point());
      sv_first->mark() = sfc->mark().boolean();
      sv_first->set_index(SE2i[sfc]);
      ++sfc;
      sv_prev = sv =
        SM.new_svertex(sface2vertex[sfc->twin()->incident_sface()]->point()-v->point());
      sv->mark() = sfc->mark().boolean();
      sv->set_index(SE2i[sfc]);
      se_first = se_prev = SM.new_shalfedge_pair(sv_first, sv_prev);
      se_first->mark() = se_first->twin()->mark() = sfc->source()->mark().boolean();
      se_first->set_index(SV2i[sfc->source()].first);
      se_first->twin()->set_index(SV2i[sfc->source()].second);
      se_first->circle() = Sphere_circle(sv_first->point(), sv->point());
      se_first->circle() = normalized(se_first->circle());
      se_first->twin()->circle() = se_first->circle().opposite();

      ++sfc;
      CGAL_For_all(sfc,sfcend) {
        sv = SM.new_svertex(sface2vertex[sfc->twin()->incident_sface()]->point()-v->point());
        sv->mark() = sfc->mark().boolean();
        sv->set_index(SE2i[sfc]);
        se = SM.new_shalfedge_pair(sv_prev, sv);
        se->mark() = se->twin()->mark() = sfc->source()->mark().boolean();
        se->set_index(SV2i[sfc->source()].first);
        se->twin()->set_index(SV2i[sfc->source()].second);
        se->circle() = Sphere_circle(sv_prev->point(), sv->point());
        se->circle() = normalized(se->circle());
        se->twin()->circle() = se->circle().opposite();
        se->sprev() = se_prev;
        se_prev->snext() = se;
        sv_prev = sv;
        se_prev = se;
      }

      se = SM.new_shalfedge_pair(sv_prev, sv_first);
      se->mark() = se->twin()->mark() = sfc->source()->mark().boolean();
      se->set_index(SV2i[sfc->source()].first);
      se->twin()->set_index(SV2i[sfc->source()].second);
      se->circle() = Sphere_circle(sv_prev->point(), sv_first->point());
      se->circle() = normalized(se->circle());
      se->twin()->circle() = se->circle().opposite();
      se->sprev() = se_prev;
      se_prev->snext() = se;
      se_first->sprev() = se;
      se->snext() = se_first;

      SFace_handle sf0 = SM.new_sface();
      SFace_handle sf1 = SM.new_sface();
      sf0->mark() = false;
      sf1->mark() = true;
      SM.link_as_face_cycle(se,sf0);
      SM.link_as_face_cycle(se->twin(),sf1);
    }
  }

  void create_single_vertex(SNC_structure& snc) {
    Vertex_handle v =
      snc.new_vertex(G.sfaces_begin()->mark().point(),
                     G.sfaces_begin()->mark().boolean());
    SM_decorator SM(&*v);
    SFace_handle sf = SM.new_sface();
    sf->mark() = false;
  }

  void create_single_edge(SNC_structure& snc) {
    SHalfloop_const_handle slc = G.shalfloop();
    Vertex_handle v0 =
      snc.new_vertex(slc->incident_sface()->mark().point(),
                     slc->incident_sface()->mark().boolean());
    SM_decorator SM0(&*v0);
    SVertex_handle sv0 = SM0.new_svertex();
    sv0->point() = slc->circle().orthogonal_vector();
    sv0->mark() = true;
    SFace_handle sf0 = SM0.new_sface();
    sf0->mark() = false;
    SM0.link_as_isolated_vertex(sv0, sf0);

    Vertex_handle v1 =
      snc.new_vertex(slc->twin()->incident_sface()->mark().point(),
                     slc->twin()->incident_sface()->mark().boolean());
    SM_decorator SM1(&*v1);
    SVertex_handle sv1 = SM1.new_svertex();
    sv1->point() = sv0->point().antipode();
    sv1->mark() = true;
    SFace_handle sf1 = SM1.new_sface();
    sf1->mark() = false;
    SM1.link_as_isolated_vertex(sv1, sf1);
  }

  void create_single_facet(SNC_structure& snc) {
       CGAL::Unique_hash_map<SHalfedge_const_handle, int> SE2i;
    SHalfedge_const_iterator sei;
    CGAL_forall_sedges(sei, G) {
      SE2i[sei] = Index_generator::get_unique_index();
      SE2i[sei->twin()] = SE2i[sei];
    }

    CGAL::Unique_hash_map
      <SVertex_const_handle, std::pair<int, int> > SV2i;
    SVertex_const_iterator svi;
    CGAL_forall_svertices(svi, G)
      SV2i[svi] = std::pair<int, int>
      (Index_generator::get_unique_index(),
       Index_generator::get_unique_index());

    CGAL::Unique_hash_map<SFace_const_handle, Vertex_handle> sface2vertex;
    SFace_const_iterator sfi;
    for(sfi = G.sfaces_begin(); sfi != G.sfaces_end(); ++sfi) {
      sface2vertex[sfi] = snc.new_vertex(sfi->mark().point(),
                                         sfi->mark().boolean());
    }

    for(sfi = G.sfaces_begin(); sfi != G.sfaces_end(); ++sfi) {
      Vertex_handle v = sface2vertex[sfi];
      SM_decorator SM(&*v);

      SHalfedge_const_handle sec = sfi->sface_cycles_begin();
      SHalfedge_around_sface_const_circulator sfc(sec);

      SVertex_handle sv0 =
        SM.new_svertex(sface2vertex[sfc->twin()->incident_sface()]->point()-v->point());
      sv0->mark() = sfc->mark().boolean();
      sv0->set_index(SE2i[sfc]);
      ++sfc;
      SVertex_handle sv1 =
        SM.new_svertex(sface2vertex[sfc->twin()->incident_sface()]->point()-v->point());
      sv1->mark() = sfc->mark().boolean();
      sv1->set_index(SE2i[sfc]);

      SHalfedge_handle se = SM.new_shalfedge_pair(sv0, sv1);
      se->mark() = se->twin()->mark() = sfc->source()->mark().boolean();
      se->set_index(SV2i[sfc->source()].first);
      se->twin()->set_index(SV2i[sfc->source()].second);
      se->circle() = Sphere_circle(sv0->point(), sv1->point());
      se->circle() = normalized(se->circle());
      se->twin()->circle() = se->circle().opposite();

      SFace_handle sf = SM.new_sface();
      sf->mark() = false;
      SM.link_as_face_cycle(se, sf);
    }
  }

  void operator()(SNC_structure& snc) {
    snc.clear();

    if(G.number_of_sfaces() == 1)
      create_single_vertex(snc);
    else if(G.number_of_sfaces() == 2)
      create_single_edge(snc);
    else if(G.number_of_svertices() == 2)
      create_single_facet(snc);
    else
      create_solid(snc);
  }


};

} //namespace CGAL
#endif // CGAL_MS3_GAUSSIAN_MAP_TO_NEF_3_H
