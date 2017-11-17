// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
