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
#ifndef CGAL_CD3_RAY_HIT_GENERATOR2_H
#define CGAL_CD3_RAY_HIT_GENERATOR2_H

#include <CGAL/license/Convex_decomposition_3.h>


#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Convex_decomposition_3/SM_walls.h>
#include <CGAL/Convex_decomposition_3/Ray_hit_generator.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 233
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template<typename Nef_>
class Ray_hit_generator2 : public Ray_hit_generator<Nef_> {
  typedef Nef_                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_and_PL    SNC_and_PL;
  typedef typename Nef_polyhedron::SNC_structure SNC_structure;
  typedef Ray_hit_generator<Nef_polyhedron>      Base;

  typedef typename SNC_structure::Sphere_map     Sphere_map;
  typedef CGAL::SM_walls<Sphere_map>             SM_walls;

  typedef typename Base::Ray_3                   Ray_3;
  typedef typename Base::Vector_3                Vector_3;
  typedef typename Base::Sphere_point            Sphere_point;
  typedef typename Base::Vertex_handle           Vertex_handle;
  typedef typename Base::SVertex_handle          SVertex_handle;
  typedef typename Base::Halfedge_handle         Halfedge_handle;

  Vertex_handle vs;
  bool edge_splitted;
  Halfedge_handle second_half;

  Vertex_handle v_new;
  bool vertex_added;

 public:
  Ray_hit_generator2(Vector_3 d, Vertex_handle v)
    : Base(d), vs(v), edge_splitted(false), vertex_added(false) {}

  void handle_splits(Halfedge_handle e, SVertex_handle svf, SVertex_handle svb) override {
      edge_splitted = true;
      if(e->source()->point() < e->twin()->source()->point())
        second_half = svf;
      else
        second_half = svb;

      CGAL_NEF_TRACEN("new edge " << e->source()->point() << "->" << e->twin()->source()->point());
      CGAL_NEF_TRACEN("new edge " << svf->source()->point() << "->" << svf->twin()->source()->point());
      CGAL_NEF_TRACEN("second_half " << second_half->source()->point()
                      << "->" << second_half->twin()->source()->point());

      vertex_added = true;
  }

  void operator()(SNC_and_PL& sncpl) override {
    Base::sncp = sncpl.sncp;
    Base::pl = sncpl.pl;

    edge_splitted = false;
    vertex_added = false;

    CGAL_NEF_TRACEN("ray hit 2: " << vs->point()
                    << " (" << dir << ")");
    SM_walls smw(&*vs);
    SVertex_handle sv1, sv2;
    if(smw.need_to_shoot(Sphere_point(Base::dir),sv1)) {
      Ray_3 r(vs->point(), Base::dir);
      v_new = Base::create_vertex_on_first_hit(r);
      SM_walls smw(&*v_new);
      sv2 = smw.add_ray_svertex(Sphere_point(-Base::dir));
      CGAL_NEF_TRACEN("sv1 " << sv1->source()->point() << "( " << sv1->point() << ")");
      CGAL_NEF_TRACEN("sv2 " << sv2->source()->point() << "( " << sv2->point() << ")");
      sv1->twin() = sv2; // TODO: why is this necessary?
      sv2->twin() = sv1; // these edges should not go into the Edge_sorter
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
      sv2->set_index(sv1->new_index());
#endif
    }
  }

  bool split_edge(Halfedge_handle& e) {
    if(edge_splitted) {
      e = second_half;
      return true;
    }
    return false;
  }

  bool added_vertex_on_edge(Vertex_handle& v) {
    if(vertex_added) {
      v = v_new;
      return true;
    }
    return false;
  }
};

} //namespace CGAL
#endif //CGAL_CD3_RAY_HIT_GENERATOR2_H
