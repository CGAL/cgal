// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
// $Id$
// 
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
#define CGAL_MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>

CGAL_BEGIN_NAMESPACE 

/// NGHK: Is the converter needed or do we just use the Cartesian_converter
template <class Triangulation_3, class HalfedgeDS, class Converter_ >
class Marching_tetrahedra_traits_skin_surface_3 {
public:
  typedef Triangulation_3                              Triangulation;
  typedef HalfedgeDS                                   Halfedge_DS;
  typedef Converter_                                   Converter;

  typedef typename Triangulation_3::Cell_handle        Cell_handle;
  

  typedef typename HalfedgeDS::Traits::Point_3         HDS_point;
  typedef typename HDS_point::R::RT                    HDS_RT;

  Marching_tetrahedra_traits_skin_surface_3() {
  }

  // These two functions are required by the marching tetrahedra algorithm
  Sign sign(Cell_handle const ch, int i) const {
    return ch->vertex(i)->sign();
  }
  HDS_point intersection(Cell_handle const ch, int i, int j) const {
    // Precondition: ch is not an infinite cell: their surface is not set
    HDS_point p1 = converter(static_cast<Triang_point>(ch->vertex(i)->point()));
    HDS_point p2 = converter(static_cast<Triang_point>(ch->vertex(j)->point()));
    return ch->surf->to_surface(p1, p2);
  }

private:
  typedef typename Triangulation::Geom_traits::Point_3 Triang_point;
  
  // Additional functions, not belonging to the traits concept:
  template <class Point>
  HDS_RT value(const Cell_handle &ch, const Point &p) const {
    return ch->surf->value(converter(p));
  }

  Converter converter;
};

CGAL_END_NAMESPACE 

#endif // CGAL_MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
