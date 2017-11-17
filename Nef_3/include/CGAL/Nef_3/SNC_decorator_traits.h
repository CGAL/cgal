// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF3_SNC_DECORATOR_TRAITS_H
#define CGAL_NEF3_SNC_DECORATOR_TRAITS_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_S2/SM_decorator_traits.h>
#include <CGAL/Nef_S2/SM_decorator.h>

namespace CGAL {

template <class Refs_>
class SNC_decorator_traits : public CGAL::SM_decorator_traits<Refs_> {
  typedef Refs_ Refs;
  typedef typename Refs::Sphere_map Sphere_map;
 public:
  typedef CGAL::SM_decorator<Sphere_map> SM_decorator;

  typedef typename Refs::Vertex_handle Vertex_handle;
  typedef typename Refs::Halfedge_handle Halfedge_handle;
  typedef typename Refs::Halffacet_handle Halffacet_handle;
  typedef typename Refs::Volume_handle Volume_handle;
  typedef typename Refs::SVertex_handle SVertex_handle;
  typedef typename Refs::SHalfedge_handle SHalfedge_handle;
  typedef typename Refs::SHalfloop_handle SHalfloop_handle;
  typedef typename Refs::SFace_handle SFace_handle;

  typedef typename Refs::Halffacet_triangle_handle 
    Halffacet_triangle_handle;

  typedef typename Refs::Vertex_iterator Vertex_iterator;
  typedef typename Refs::Halfedge_iterator Halfedge_iterator;
  typedef typename Refs::Halffacet_iterator Halffacet_iterator; 
  typedef typename Refs::Volume_iterator Volume_iterator;
  typedef typename Refs::SVertex_iterator SVertex_iterator;  
  typedef typename Refs::SHalfedge_iterator SHalfedge_iterator;
  typedef typename Refs::SHalfloop_iterator SHalfloop_iterator;
  typedef typename Refs::SFace_iterator SFace_iterator;  

  typedef typename Refs::SHalfedge_around_svertex_circulator 
    SHalfedge_around_svertex_circulator;
  typedef typename Refs::SHalfedge_around_sface_circulator 
    SHalfedge_around_sface_circulator;
  typedef typename Refs::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename Refs::SHalfedge_around_facet_circulator 
    SHalfedge_around_facet_circulator;
  typedef typename Refs::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename Refs::Shell_entry_iterator Shell_entry_iterator;
};

template <class Refs_>
class SNC_decorator_const_traits {
  typedef Refs_ Refs;
  typedef typename Refs::Sphere_map Sphere_map;
 public:
  typedef CGAL::SM_const_decorator<Sphere_map> SM_decorator;

  typedef typename Refs::Vertex_const_handle Vertex_handle;
  typedef typename Refs::Halfedge_const_handle Halfedge_handle;
  typedef typename Refs::Halffacet_const_handle Halffacet_handle;
  typedef typename Refs::Volume_const_handle Volume_handle;
  typedef typename Refs::SVertex_const_handle SVertex_handle;
  typedef typename Refs::SHalfedge_const_handle SHalfedge_handle;
  typedef typename Refs::SHalfloop_const_handle SHalfloop_handle;
  typedef typename Refs::SFace_const_handle SFace_handle;

  typedef typename Refs::Halffacet_triangle_const_handle 
    Halffacet_triangle_handle;

  typedef typename Refs::Vertex_const_iterator Vertex_iterator;
  typedef typename Refs::Halfedge_const_iterator Halfedge_iterator;
  typedef typename Refs::Halffacet_const_iterator Halffacet_iterator; 
  typedef typename Refs::Volume_const_iterator Volume_iterator;
  typedef typename Refs::SVertex_const_iterator SVertex_iterator;  
  typedef typename Refs::SHalfedge_const_iterator SHalfedge_iterator;
  typedef typename Refs::SHalfloop_const_iterator SHalfloop_iterator;
  typedef typename Refs::SFace_const_iterator SFace_iterator;  

  typedef typename Refs::SHalfedge_around_svertex_const_circulator 
    SHalfedge_around_svertex_circulator;
  typedef typename Refs::SHalfedge_around_sface_const_circulator 
    SHalfedge_around_sface_circulator;
  typedef typename Refs::SFace_cycle_const_iterator SFace_cycle_iterator;
  typedef typename Refs::SHalfedge_around_facet_const_circulator 
    SHalfedge_around_facet_circulator;
  typedef typename Refs::Halffacet_cycle_const_iterator Halffacet_cycle_iterator;
  typedef typename Refs::Shell_entry_const_iterator Shell_entry_iterator;
};

} //namespace CGAL
#endif // CGAL_NEF3_SNC_DECORATOR_TRAITS_H
