// Copyright (c) 2005 RuG (Netherlands)
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
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
#define MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>

CGAL_BEGIN_NAMESPACE 

template <class Triangulation_3, class HalfedgeDS, class Converter_ >
class Marching_tetrahedra_traits_skin_surface_3 {
public:
  typedef Triangulation_3                              Triangulation;
  typedef HalfedgeDS                                   Halfedge_DS;
  typedef Converter_                                   Converter;

  typedef typename Triangulation::Vertex_handle        Vertex_handle;
  typedef typename Triangulation::Edge                 Edge;
  typedef typename Triangulation_3::Cell_handle        Cell_handle;
  typedef typename Triangulation::Geom_traits::Point_3 Triang_point;

  typedef typename HalfedgeDS::Traits                  HDS_K;
  typedef typename HDS_K::Point_3                      HDS_point;
  typedef typename HDS_point::R::RT                    HDS_rt;

  Marching_tetrahedra_traits_skin_surface_3(HDS_rt iso_value=0)
    : iso_value(iso_value) {
  }

  // These two functions are required by the marching tetrahedra algorithm
  Sign sign(Cell_handle const ch, int i) const {
    return CGAL_NTS sign(value(ch,ch->vertex(i)->point()) - iso_value);
  }
  HDS_point intersection(Cell_handle const ch, int i, int j) const {
    // Precondition: ch is not an infinite cell: their surface is not set
    HDS_point p1 = converter(static_cast<Triang_point>(ch->vertex(i)->point()));
    HDS_point p2 = converter(static_cast<Triang_point>(ch->vertex(j)->point()));
    return ch->surf->to_surface(p1, p2);
  }

private:
  // Additional functions, not belonging to the traits concept:
  HDS_rt value(const Cell_handle &ch, const HDS_point &p) const {
    return ch->surf->value(p);
  }
  HDS_rt value(const Cell_handle &ch, const Triang_point &p) const {
    return ch->surf->value(converter(p));
  }

  Converter converter;
  HDS_rt iso_value;
};

CGAL_END_NAMESPACE 

#endif // MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
