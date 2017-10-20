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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_CD3_YVERTICAL_WALL_BUILDER_H
#define CGAL_CD3_YVERTICAL_WALL_BUILDER_H

#include <CGAL/license/Convex_decomposition_3.h>


#include<CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Convex_decomposition_3/Single_wall_creator3.h>
#include <CGAL/Convex_decomposition_3/External_structure_builder.h>

namespace CGAL {

template<typename Nef_>
class YVertical_wall_builder : public Modifier_base<typename Nef_::SNC_and_PL> {

  typedef Nef_                                            Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_and_PL             SNC_and_PL;
  typedef typename Nef_polyhedron::SNC_structure          SNC_structure;
  typedef CGAL::SNC_decorator<SNC_structure>              SNC_decorator;
  typedef CGAL::SNC_decorator<SNC_structure>     Base;
  typedef typename CGAL::Single_wall_creator3<Nef_polyhedron> Single_wall3;
  typedef typename CGAL::External_structure_builder<Nef_polyhedron> External_structure_builder;

  typedef typename SNC_structure::Vertex_handle           Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle         Halfedge_handle;
  typedef typename SNC_structure::Halffacet_handle        Halffacet_handle;
  typedef typename SNC_structure::SHalfedge_handle        SHalfedge_handle;
  typedef typename SNC_structure::SHalfloop_handle        SHalfloop_handle;
  typedef typename SNC_structure::SFace_handle            SFace_handle;

  typedef typename SNC_structure::Volume_iterator         Volume_iterator;
  typedef typename SNC_structure::Halfedge_iterator       Halfedge_iterator;
  typedef typename SNC_structure::SFace_iterator       SFace_iterator;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator
    SHalfedge_around_svertex_circulator;

  typedef typename SNC_structure::Vector_3                Vector_3;
  typedef typename SNC_structure::Point_3                 Point_3;

  typedef typename SNC_structure::Sphere_point            Sphere_point;
  typedef typename SNC_structure::Sphere_circle           Sphere_circle;
  typedef typename SNC_structure::Sphere_segment          Sphere_segment;

  typedef typename std::list<Halfedge_handle>             Edge_list;
 public:
  typedef typename std::list<Halfedge_handle>::iterator   Vertical_redge_iterator;

  Edge_list redges;

 public:
  YVertical_wall_builder() {}
    
  void operator()(SNC_and_PL& sncpl) {
    SNC_structure* sncp(sncpl.sncp);

    SFace_iterator sfi;
    CGAL_forall_sfaces(sfi, *sncp)
      if(sncp->is_boundary_object(sfi))
	sncp->undef_boundary_item(sfi);

    Halfedge_iterator ei;
    CGAL_forall_halfedges(ei, *sncp) {
      if(ei->point() != Sphere_point(1,0,0)) continue;
      SHalfedge_around_svertex_circulator 
	svc(ei->out_sedge()), send(svc);
      CGAL_For_all(svc, send) {
	if(!svc->incident_sface()->mark()) continue;
	if(!CGAL::is_reflex_sedge_in_any_direction<SNC_structure>(svc))
	  continue;
	redges.push_back(ei);
	break;
      }
    }
    
    Vertical_redge_iterator vri;
    for(vri = redges_begin(); vri != redges_end(); ++vri) {
      Halfedge_handle ei(*vri);
      SHalfedge_around_svertex_circulator 
	svc(ei->out_sedge()), send(svc);
      CGAL_For_all(svc, send) {
	if(!svc->incident_sface()->mark()) continue;
	if(!CGAL::is_reflex_sedge_in_any_direction<SNC_structure>(svc))
	  continue;
	Single_wall3 W(svc);
	W(sncpl);
	break;
      }
    }

  }

  Vertical_redge_iterator redges_begin() { return redges.begin(); }
  Vertical_redge_iterator redges_end()   { return redges.end(); }
};

} //namespace CGAL
#endif // CGAL_CD3_YVERTICAL_WALL_BUILDER_H
