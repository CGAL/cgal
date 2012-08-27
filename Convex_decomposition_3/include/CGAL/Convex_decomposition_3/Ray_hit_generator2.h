// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_CD3_RAY_HIT_GENERATOR2_H
#define CGAL_CD3_RAY_HIT_GENERATOR2_H

#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_intersection.h>
#include <CGAL/Convex_decomposition_3/SM_walls.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 233
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template<typename Nef_>
class Ray_hit_generator2 : public Modifier_base<typename Nef_::SNC_and_PL> {
  
  typedef Nef_                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_and_PL    SNC_and_PL;
  typedef typename Nef_polyhedron::SNC_structure SNC_structure;
  typedef typename Nef_polyhedron::Items         Items;
  typedef CGAL::SNC_decorator<SNC_structure>     Base;
  typedef CGAL::SNC_point_locator<Base>          SNC_point_locator;
  typedef CGAL::SNC_intersection<SNC_structure>  SNC_intersection;
  typedef CGAL::SNC_constructor<Items, SNC_structure>   
    SNC_constructor;

  typedef typename SNC_structure::Sphere_map     Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>         SM_decorator;  
  typedef CGAL::SM_point_locator<SM_decorator>   SM_point_locator; 
  typedef CGAL::SM_walls<Sphere_map>             SM_walls;

  typedef typename Base::Segment_3               Segment_3;
  typedef typename Base::Point_3                 Point_3;
  typedef typename Base::Ray_3                   Ray_3;
  typedef typename Base::Vector_3                Vector_3;
  typedef typename Base::Sphere_point            Sphere_point;
  typedef typename Base::Vertex_handle           Vertex_handle;
  typedef typename Base::SVertex_handle           SVertex_handle;
  typedef typename Base::Halfedge_handle         Halfedge_handle;
  typedef typename Base::Halffacet_handle        Halffacet_handle;
  typedef typename Base::Object_handle           Object_handle;

  typedef typename Base::Vertex_iterator         Vertex_iterator;
  typedef typename Base::SVertex_iterator         SVertex_iterator;

  Vector_3 dir;
  Vertex_handle vs;
  SNC_structure* sncp;
  SNC_point_locator* pl;

  bool edge_splitted;
  Halfedge_handle second_half;

  Vertex_handle v_new;
  bool vertex_added;

 public:
  Ray_hit_generator2(Vector_3 d, Vertex_handle v) : dir(d), vs(v) {}
      
  Vertex_handle create_vertex_on_first_hit(const Ray_3& r) {

    CGAL_NEF_TRACEN("shoot ray in SNC " << r);

    //    CGAL_NEF_SETDTHREAD(503*509);
    Object_handle o = pl->shoot(r);
    //    CGAL_NEF_SETDTHREAD(1);

    Vertex_handle v;
    if(assign(v, o)) {
      CGAL_NEF_TRACEN("Found vertex " << v->point());
      return v;
    }

    Point_3 ip;
    SNC_intersection I;
    SNC_constructor C(*sncp);

    Halfedge_handle e;
    if(assign(e, o)) {
       CGAL_NEF_TRACEN("Found edge " << e->source()->point() 
		       << "->" << e->twin()->source()->point());
      Segment_3 seg(e->source()->point(), e->twin()->source()->point());
      I.does_intersect_internally(r, seg, ip);
      ip = normalized(ip);
      v = C.create_from_edge(e,ip);
      pl->add_vertex(v);

      CGAL_NEF_TRACEN("new vertex " << ip);

      SVertex_iterator svi = v->svertices_begin();
      SVertex_handle svf = svi;
      SVertex_handle svb = ++svi;

      if(svf->point() == e->point()) {
	svb->twin() = e;
	svf->twin() = e->twin();
	e->twin()->twin() = svf;
	e->twin() = svb;
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
	svb->set_index(e->get_index());
	svf->set_index();
	svf->twin()->set_index(svf->get_index());
#endif
      } else {
	svf->twin() = e;
	svb->twin() = e->twin();
	e->twin()->twin() = svb;
	e->twin() = svf;
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
	svf->set_index(e->get_index());
	svb->set_index();
	svb->twin()->set_index(svb->get_index());
#endif
      }

      pl->add_edge(svf);
      pl->add_edge(svb);

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
      return v;
    }

    Halffacet_handle f;
    if(assign(f, o)) {
      CGAL_NEF_TRACEN("Found facet ");
      I.does_intersect_internally(r, f, ip);
      ip = normalized(ip);
      v = C.create_from_facet(f,ip);
      pl->add_vertex(v);
      CGAL_NEF_TRACEN("new vertex " << ip);

      return v;
    }

    CGAL_error_msg( "ray should hit vertex, edge, or facet");
    return Vertex_handle();
  }

  void operator()(SNC_and_PL& sncpl) {
    sncp = sncpl.sncp;
    pl = sncpl.pl;

    edge_splitted = false;
    vertex_added = false;
    
    CGAL_NEF_TRACEN("ray hit 2: " << vs->point()
		    << " (" << dir << ")");
    SM_walls smw(&*vs);
    SVertex_handle sv1, sv2;
    if(smw.need_to_shoot(Sphere_point(dir),sv1)) {
      Ray_3 r(vs->point(), dir);
      v_new = create_vertex_on_first_hit(r);
      SM_walls smw(&*v_new);
      sv2 = smw.add_ray_svertex(Sphere_point(-dir));
      CGAL_NEF_TRACEN("sv1 " << sv1->source()->point() << "( " << sv1->point() << ")");
      CGAL_NEF_TRACEN("sv2 " << sv2->source()->point() << "( " << sv2->point() << ")");
      sv1->twin() = sv2; // TODO: why is this necessary?
      sv2->twin() = sv1; // these edges should not go into the Edge_sorter
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
      sv1->set_index();
      sv2->set_index(sv1->get_index());
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
