// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
#define CGAL_MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H

#include <CGAL/license/Skin_surface_3.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Cartesian_converter.h>

namespace CGAL {

/// NGHK: Is the converter needed or do we just use the Cartesian_converter
template <class SkinSurface_3,
          class Vertex_iterator,
          class Cell_iterator,
          class HalfedgeDS>
class Marching_tetrahedra_traits_skin_surface_3
{
public:
  typedef HalfedgeDS                                   Halfedge_DS;

  typedef typename HalfedgeDS::Traits::Point_3         HDS_point;
  typedef typename HDS_point::R::RT                    HDS_RT;

  typedef typename SkinSurface_3::Bare_point           Skin_point;

  Marching_tetrahedra_traits_skin_surface_3(const SkinSurface_3 &ss_3)
    : ss_3(ss_3)
  { }

  // These two functions are required by the marching tetrahedra algorithm
  Sign sign(const Cell_iterator ch, int i) const {
    return ss_3.sign(ch->vertex(i));
  }

  HDS_point intersection(Cell_iterator const ch, int i, int j) const
  {
    // Precondition: ch is not an infinite cell: their surface is not set
    Skin_point p;
    ss_3.intersect(ch, i, j, p);

    return
      Cartesian_converter<typename Skin_point::R, typename HDS_point::R>()(p);
  }

  const SkinSurface_3 &ss_3;
};

} //namespace CGAL

#endif // CGAL_MARCHING_TETRAHEDRA_TRAITS_SKIN_SURFACE_3_H
