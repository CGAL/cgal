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

#ifndef CGAL_MAKE_SKIN_SURFACE_MESH_3_H
#define CGAL_MAKE_SKIN_SURFACE_MESH_3_H

#include <CGAL/license/Skin_surface_3.h>


#include <CGAL/Skin_surface_traits_3.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>

#include <CGAL/make_union_of_balls_3.h>

namespace CGAL {

template <class WP_iterator,
          class Polyhedron_3>
void make_skin_surface_mesh_3(Polyhedron_3 &p,
                              WP_iterator begin, WP_iterator end,
                              double shrink_factor=.5,
                              int nSubdivisions=0,
                              bool grow_balls=true)
{
  if (shrink_factor == 1) {
    make_union_of_balls_mesh_3(p,begin,end,nSubdivisions);
  }

  typedef typename WP_iterator::value_type               Weighted_point;
  typedef typename Kernel_traits<Weighted_point>::Kernel K;

  typedef Skin_surface_traits_3<K>                       Traits;
  typedef Skin_surface_3<Traits>                         Skin_surface;

  Skin_surface skin_surface(begin, end, shrink_factor, grow_balls);

  mesh_skin_surface_3(skin_surface, p);

  subdivide_skin_surface_mesh_3(skin_surface, p, nSubdivisions);
}



} //namespace CGAL

#endif // CGAL_MAKE_SKIN_SURFACE_MESH_3_H
