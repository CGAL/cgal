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
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
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
#ifndef CGAL_NEF_SHELL_TO_NEF_3_H
#define CGAL_NEF_SHELL_TO_NEF_3_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Circulator_project.h>
#include <CGAL/normal_vector_newell_3.h>
#include <CGAL/Nef_S2/SM_point_locator.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 29
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template <class SNC_structure>
class Shell_to_nef_3
{
  typedef typename SNC_structure::SM_decorator       SM_decorator;
  typedef typename SNC_structure::Vertex_handle      Vertex_handle;
  typedef typename SNC_structure::SVertex_handle     SVertex_handle;
  typedef typename SNC_structure::SHalfedge_handle   SHalfedge_handle;
  typedef typename SNC_structure::SFace_handle       SFace_handle;
  typedef typename SNC_structure::Point_3            Point_3;
  typedef typename SNC_structure::Plane_3            Plane_3;
  typedef typename SNC_structure::Sphere_point       Sphere_point;
  typedef typename SNC_structure::Sphere_segment     Sphere_segment;
  typedef typename SNC_structure::Sphere_circle      Sphere_circle;

  typedef typename SNC_structure::Vertex_const_handle
    Vertex_const_handle;
  typedef typename SNC_structure::Halfedge_const_handle
    Halfedge_const_handle;
  typedef typename SNC_structure::Halffacet_const_handle
    Halffacet_const_handle;
  typedef typename SNC_structure::SHalfedge_const_handle
    SHalfedge_const_handle;
  typedef typename SNC_structure::SHalfloop_const_handle
    SHalfloop_const_handle;
  typedef typename SNC_structure::SFace_const_handle 
    SFace_const_handle;
  typedef typename SNC_structure::SFace_cycle_const_iterator
    SFace_cycle_const_iterator;
  typedef typename SNC_structure::SHalfedge_around_sface_const_circulator
    SHalfedge_around_sface_const_circulator;

 private:
  SNC_structure& S;
 public:

  Shell_to_nef_3(SNC_structure& S_) : S(S_) {}
  
  void visit(Vertex_const_handle ) {}
  void visit(Halfedge_const_handle ) {}
  void visit(Halffacet_const_handle ) {}
  void visit(SHalfedge_const_handle ) {}
  void visit(SHalfloop_const_handle ) {}

  void visit(SFace_const_handle sf) {

    SFace_cycle_const_iterator sfci =
      sf->sface_cycles_begin();
    if(sfci.is_shalfloop())
      return;
    CGAL_assertion(sfci.is_shalfedge());

    Vertex_const_handle pv = sf->center_vertex();
    Vertex_handle nv = S.new_vertex();
    nv->point() = pv->point();
    nv->mark() = true;
      
    SM_decorator SM(&*nv);
    SHalfedge_around_sface_const_circulator 
      pe(sfci), pe_prev(pe), pend(pe);
      
    SVertex_handle sv_0 = 
      SM.new_svertex(pe->source()->point());
    sv_0->mark() = true;
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
    sv_0->set_index(pe->source()->get_index());
#endif
    ++pe;  
    SVertex_handle sv_prev = sv_0;
      
    do {
      SVertex_handle sv = SM.new_svertex(pe->source()->point());
      sv->mark() = true;
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
      sv->set_index(pe->source()->get_index());
#endif
      
      SHalfedge_handle e = SM.new_shalfedge_pair(sv_prev, sv);
      e->circle() = pe_prev->circle();
      e->twin()->circle() = pe_prev->twin()->circle();
      e->mark() = e->twin()->mark() = true;
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
      e->set_index(pe_prev->get_index());
      e->twin()->set_index(pe_prev->twin()->get_index());
#endif	
      sv_prev = sv;
      pe_prev = pe;
      ++pe;
    }
    while( pe != pend );
      
    SHalfedge_handle e;
    e = SM.new_shalfedge_pair(sv_prev, sv_0);
    e->circle() = pe_prev->circle();
    e->twin()->circle() = pe_prev->twin()->circle();
    e->mark() = e->twin()->mark() = true;
#ifndef CGAL_NEF_NO_INDEXED_ITEMS
    e->set_index(pe_prev->get_index());
    e->twin()->set_index(pe_prev->twin()->get_index());    
#endif
    // create faces
    SFace_handle fext = SM.new_sface();
    SM.link_as_face_cycle(e->twin(), fext);
    fext->mark() = false;
    
    SFace_handle fint = SM.new_sface();
    SM.link_as_face_cycle(e, fint);
    fint->mark() = false;
    
    SM.check_integrity_and_topological_planarity();   
  }
};

template <class Nef_polyhedron>
void shell_to_nef_3(const Nef_polyhedron& N,
		    typename Nef_polyhedron::SFace_const_handle sf,
		    typename Nef_polyhedron::SNC_structure& S) 
{
  Shell_to_nef_3<typename Nef_polyhedron::SNC_structure> s2n(S);
  N.visit_shell_objects(sf, s2n);
}

} //namespace CGAL

#endif //CGAL_NEF_SHELL_TO_NEF_3_H
