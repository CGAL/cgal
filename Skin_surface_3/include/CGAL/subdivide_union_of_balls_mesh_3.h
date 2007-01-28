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
// $URL: svn+ssh://nicokruithof@scm.gforge.inria.fr/svn/cgal/trunk/Skin_surface_3/include/CGAL/subdivide_union_of_balls_mesh_3.h $
// $Id: subdivide_union_of_balls_mesh_3.h 35761 2007-01-19 16:53:56Z nicokruithof $
// 
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_SUBDIVIDE_UNION_OF_BALLS_MESH_3_H
#define CGAL_SUBDIVIDE_UNION_OF_BALLS_MESH_3_H

#include <CGAL/Skin_surface_refinement_policy_3.h>
#include <CGAL/Polyhedron_3.h>

CGAL_BEGIN_NAMESPACE

template <class UnionOfBalls_3, class Polyhedron_3>
void subdivide_union_of_balls_mesh_3(
          const UnionOfBalls_3 &skin,
          Polyhedron_3 &p, 
          int nSubdiv = 1) {
  while (nSubdiv > 0) {
    skin.subdivide_mesh_3(p);
    nSubdiv--;
  }
}

CGAL_END_NAMESPACE

#endif // CGAL_SUBDIVIDE_UNION_OF_BALLS_MESH_3_H
