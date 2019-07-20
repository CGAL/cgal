// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// WARRANTY OF DESISGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
#ifndef CGAL_COMBINE_WITH_HALFSPACE_H
#define CGAL_COMBINE_WITH_HALFSPACE_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_S2/Normalizing.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_3/SNC_const_decorator.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_3/SNC_SM_overlayer.h>
#include <CGAL/Nef_3/SNC_constructor.h>
#include <CGAL/Nef_3/SNC_external_structure.h>
#include <CGAL/Nef_3/ID_support_handler.h>
#include <CGAL/Nef_3/Binary_operation.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 19
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <typename Map, typename SNC_point_locator>
class Combine_with_halfspace : public SNC_decorator<Map> { 
 public:
  typedef Map SNC_structure;
  typedef typename SNC_structure::Items                Items;
  typedef typename Map::Sphere_map                     Sphere_map;
  typedef CGAL::Combine_with_halfspace
    <SNC_structure, SNC_point_locator>  Self;
  typedef CGAL::SNC_decorator<SNC_structure>           Base;
  typedef Base                                         SNC_decorator;
  typedef CGAL::SNC_constructor<Items, SNC_structure>  SNC_constructor;
  typedef CGAL::SNC_external_structure<Items, SNC_structure> 
    SNC_external_structure;
  typedef CGAL::SM_decorator<Sphere_map>               SM_decorator;
  typedef CGAL::SM_const_decorator<Sphere_map>         SM_const_decorator;
  typedef CGAL::Binary_operation<SNC_structure>        Binary_operation;

  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;

  typedef typename SNC_structure::Vertex_const_iterator Vertex_const_iterator;
  typedef typename SNC_structure::Halfedge_const_iterator Halfedge_const_iterator;

  typedef typename Base::SHalfedge_iterator SHalfedge_iterator;
  typedef typename Base::SHalfedge_const_iterator SHalfedge_const_iterator;
  typedef typename Base::SHalfloop_const_iterator SHalfloop_const_iterator;

  typedef typename Base::Point_3 Point_3;
  typedef typename Base::Segment_3 Segment_3;
  typedef typename Base::Plane_3 Plane_3;

  typedef typename Base::Mark Mark;

  typedef CGAL::ID_support_handler<Items, Base> Association;

  SNC_point_locator* pl;

 public:
  enum Intersection_mode { CLOSED_HALFSPACE=0, OPEN_HALFSPACE=1, PLANE_ONLY=2 };

  Combine_with_halfspace(SNC_structure& W, SNC_point_locator* pl_) 
    : Base(W), pl(pl_) {}
    
  template <typename Selection, typename Intersection_mode>
  void combine_with_halfspace(const SNC_structure& snc,
			      const Plane_3& plane,
			      const Selection& BOP,
			      Intersection_mode im) {
    
    Association A;
    SHalfedge_const_iterator sei;
    CGAL_forall_shalfedges(sei, snc)
      A.initialize_hash(sei);
    SHalfloop_const_iterator sli;
    CGAL_forall_shalfloops(sli, snc)
      A.initialize_hash(sli);
    
    int index0(Index_generator::get_unique_index());
    int index1(Index_generator::get_unique_index());
    Halffacet_handle dummy_facet = 
      this->sncp()->new_halffacet_pair(plane, im != OPEN_HALFSPACE);
	A.initialize_hash(index0);
	A.initialize_hash(index1);

    Binary_operation bo(*this->sncp());
      Vertex_const_iterator v0;
      CGAL_forall_vertices( v0, snc) {
	Oriented_side os = plane.oriented_side(v0->point());
	if(os == ON_ORIENTED_BOUNDARY) {
	  SNC_constructor C(*this->sncp());
	  Vertex_handle vp =
	    C.create_from_plane(plane, v0->point(),
				im != OPEN_HALFSPACE, 
				im != PLANE_ONLY, false);
      vp->shalfloop()->set_index_facet(dummy_facet);
      vp->shalfloop()->twin()->set_index_facet(dummy_facet->twin());
      vp->shalfloop()->set_index(index0);
      vp->shalfloop()->twin()->set_index(index1);
	  Vertex_handle vr = 
	    bo.binop_local_views(v0, vp, BOP, *this->sncp(), A);
	  this->sncp()->delete_vertex(vp);
	} else if(os == ON_NEGATIVE_SIDE && im != PLANE_ONLY) {
	  SNC_constructor C(*this->sncp());
	  Vertex_handle v1 = C.clone_SM(v0);	
	}
      }

    Halfedge_const_iterator e0;
    CGAL_forall_edges(e0, snc) {
      Segment_3 seg(e0->source()->point(), 
		    e0->twin()->source()->point());
      Object o = intersection(plane, seg);
      Point_3 ip;
      if(!assign(ip,o)) continue;
      // TODO: optimize for filtering
      ip = normalized(ip);
      if(ip == e0->source()->point() ||
	 ip == e0->twin()->source()->point()) continue;
      SNC_constructor C(*this->sncp());
      Vertex_handle vp = 
	C.create_from_plane(plane, ip,
			    im != OPEN_HALFSPACE, 
			    im != PLANE_ONLY, false);

      vp->shalfloop()->set_index_facet(dummy_facet);
      vp->shalfloop()->twin()->set_index_facet(dummy_facet->twin());
      vp->shalfloop()->set_index(index0);
      vp->shalfloop()->twin()->set_index(index1);

      Vertex_handle ve = C.create_from_edge(e0, ip);
      
      Vertex_handle vr = 
	bo.binop_local_views(ve, vp, BOP, *this->sncp(), A);
      this->sncp()->delete_vertex(vp);
      this->sncp()->delete_vertex(ve);
    }
    
    this->sncp()->delete_halffacet_pair(dummy_facet);

    SHalfedge_iterator se;
    CGAL_forall_sedges(se, *this->sncp()) {
      se->circle() = normalized(se->circle());
      se->twin()->circle() = se->circle().opposite();
    }

    SNC_external_structure es(*this->sncp(), pl);
    es.build_after_binary_operation(A);	
  }
};

} //namespace CGAL
#endif //CGAL_COMBINE_WITH_HALFSPACE_H
