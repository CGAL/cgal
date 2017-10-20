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

#ifndef CGAL_CD3_EXTERNAL_STRUCTURE_BUILDER_H
#define CGAL_CD3_EXTERNAL_STRUCTURE_BUILDER_H

#include <CGAL/license/Convex_decomposition_3.h>


#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_intersection.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 43
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template<typename Nef_>
class External_structure_builder : public Modifier_base<typename Nef_::SNC_and_PL> {
  
  typedef Nef_                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_and_PL    SNC_and_PL;
  typedef typename Nef_polyhedron::SNC_structure SNC_structure;
  typedef typename SNC_structure::Items          Items;
  typedef CGAL::SNC_decorator<SNC_structure>     Base;
  typedef CGAL::SNC_point_locator<Base>          SNC_point_locator;
  typedef CGAL::SNC_intersection<SNC_structure>  SNC_intersection;
  typedef CGAL::SNC_external_structure<Items, SNC_structure>
    SNC_external_structure;

  typedef typename SNC_structure::Sphere_map     Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>         SM_decorator;  
  typedef CGAL::SM_point_locator<SM_decorator>   SM_point_locator; 

  typedef typename Base::Segment_3               Segment_3;
  typedef typename Base::Point_3                 Point_3;
  typedef typename Base::Ray_3                   Ray_3;
  typedef typename Base::Vector_3                Vector_3;
  typedef typename Base::Sphere_point            Sphere_point;
  typedef typename Base::Sphere_circle           Sphere_circle;
  typedef typename Base::Sphere_segment          Sphere_segment;
  typedef typename Base::Vertex_handle           Vertex_handle;
  typedef typename Base::Halfedge_handle         Halfedge_handle;
  typedef typename Base::Halffacet_handle        Halffacet_handle;
  typedef typename Base::SVertex_handle          SVertex_handle;
  typedef typename Base::SHalfedge_handle        SHalfedge_handle;
  typedef typename Base::SHalfloop_handle        SHalfloop_handle;
  typedef typename Base::SFace_handle            SFace_handle;
  typedef typename Base::Object_handle           Object_handle;

  typedef typename Base::SFace_iterator          SFace_iterator;
  typedef typename Base::SHalfedge_iterator      SHalfedge_iterator;
  typedef typename Base::SFace_cycle_iterator    SFace_cycle_iterator;
  typedef typename Base::SHalfedge_around_sface_circulator
    SHalfedge_around_sface_circulator;

  Halfedge_handle ein;
  Vector_3 dir;
  
 public:
  External_structure_builder() {}

  void operator()(SNC_and_PL& sncpl) {
    //    CGAL_NEF_TRACEN(43);

    SNC_structure* sncp(sncpl.sncp);
    SNC_point_locator* pl(sncpl.pl);



    Unique_hash_map<SHalfedge_handle, SFace_handle> sedge2sface;
    /*    
    SFace_iterator sfi;
    CGAL_forall_sfaces(sfi, *sncp) {
      SFace_cycle_iterator sfc;
      for(sfc = sfi->sface_cycles_begin(); sfc != sfi->sface_cycles_end(); ++sfc) {
	if(sfc.is_shalfedge()){
	  SHalfedge_around_sface_circulator eaf(sfc), end(eaf);
	  CGAL_For_all(eaf,end) {
	    SHalfedge_handle se(eaf);
	    sedge2sface[eaf] = sfi;
	  }
	}
      }
    }

    //    CGAL::SNC_io_parser<SNC_structure> O0(std::cerr, *sncp, false);
    //    O0.print();

    SHalfedge_iterator sei;
    CGAL_forall_shalfedges(sei, *sncp) {
      SHalfedge_handle se(sei);
      if(sedge2sface[se] == SFace_handle()) {
	SM_decorator SD(&*sei->source()->source());
	SFace_handle sf_new = SD.new_sface();
	sf_new->mark() = sei->incident_sface()->mark();
	
	CGAL_NEF_TRACEN("new entry sedge " << sei->source()->point() 
			<< "->" << sei->twin()->source()->point() 
			<< " at " << sei->source()->source()->point());

	SD.link_as_face_cycle(sei, sf_new);
	
	SHalfedge_around_sface_circulator eaf(se), end(eaf);
	CGAL_For_all(eaf,end) {
	  SHalfedge_handle se(eaf);
	  sedge2sface[eaf] = sf_new;
	}	

	// TODO: relink inner sface cycles
      }
    }
    */
    SNC_point_locator* old_pl = pl;
    pl = pl->clone();
    sncpl.pl = pl;
    delete old_pl;
    SNC_external_structure C(*sncp,pl);
    C.clear_external_structure();
    C.build_external_structure();

    //    CGAL::SNC_io_parser<SNC_structure> Ox(std::cerr, *sncp, false);
    //    Ox.print();
  }
};

} //namespace CGAL
#endif //CGAL_CD3_EXTERNAL_STRUCTURE_BUILDER_H
