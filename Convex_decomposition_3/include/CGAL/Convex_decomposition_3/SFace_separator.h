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
#ifndef CGAL_CD3_SFACE_SEPARATOR_H
#define CGAL_CD3_SFACE_SEPARATOR_H

#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Convex_decomposition_3/SM_walls.h>

CGAL_BEGIN_NAMESPACE

template<typename Nef_>
class SFace_separator : public Modifier_base<typename Nef_::SNC_structure> {
  
  typedef Nef_                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_structure SNC_structure;
  typedef typename SNC_structure::Items          Items;
  typedef CGAL::SNC_decorator<SNC_structure>     Base;
  typedef CGAL::SNC_point_locator<Base>          SNC_point_locator;
  typedef CGAL::SNC_constructor<Items, SNC_structure>   
    SNC_constructor;

  typedef typename SNC_structure::Sphere_map     Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>         SM_decorator;  
  typedef CGAL::SM_point_locator<SM_decorator>   SM_point_locator; 
  typedef CGAL::SM_walls<Sphere_map>             SM_walls;

  typedef typename Base::SHalfedge_handle        SHalfedge_handle;
  typedef typename Base::SHalfloop_handle        SHalfloop_handle;
  typedef typename Base::SFace_handle            SFace_handle;

  typedef typename Base::SFace_iterator          SFace_iterator;
  typedef typename Base::SFace_cycle_iterator    SFace_cycle_iterator;

 public:
  SFace_separator() {}

  void operator()(SNC_structure& snc) {

    SFace_iterator sf;
    CGAL_forall_sfaces(sf, snc) {
      if(!sf->mark() ||
	 sf->sface_cycles_begin() == 
	 sf->sface_cycles_end()) continue;

      SM_decorator SD(&*sf->center_vertex());

      SFace_cycle_iterator sfci
	(++sf->sface_cycles_begin());
      while(sfci != sf->sface_cycles_end()) {
	SFace_handle sf_new = SD.new_sface();
	sf_new->mark() = sf->mark();
	sf_new->volume() = sf->volume();
	if(sfci.is_shalfedge()) {
	  SHalfedge_handle se = sfci;
	  SD.unlink_as_face_cycle(se);
	  SD.link_as_face_cycle(se,sf_new);
	} else if(sfci.is_shalfloop()) {
       	  SHalfloop_handle sl = sfci;
	  SD.unlink_as_loop(sl);
	  SD.link_as_loop(sl, sf_new);
	} else
	  CGAL_error_msg("there should be no isolated edges");
	sfci = ++sf->sface_cycles_begin();
      }
    }
  }
};

CGAL_END_NAMESPACE
#endif //CGAL_CD3_SFACE_SEPARATOR_H
