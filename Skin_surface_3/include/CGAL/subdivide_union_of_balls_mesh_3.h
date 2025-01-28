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

#ifndef CGAL_SUBDIVIDE_UNION_OF_BALLS_MESH_3_H
#define CGAL_SUBDIVIDE_UNION_OF_BALLS_MESH_3_H

#include <CGAL/license/Skin_surface_3.h>

#include <CGAL/Skin_surface_refinement_policy_3.h>
#include <CGAL/Polyhedron_3.h>

namespace CGAL {

template <class UnionOfBalls_3, class Polyhedron_3>
void subdivide_union_of_balls_mesh_3(const UnionOfBalls_3 &skin,
                                     Polyhedron_3 &p,
                                     int nSubdiv = 1)
{
  while (nSubdiv > 0) {
    skin.subdivide_mesh_3(p);
    nSubdiv--;
  }
}

} //namespace CGAL

#endif // CGAL_SUBDIVIDE_UNION_OF_BALLS_MESH_3_H
