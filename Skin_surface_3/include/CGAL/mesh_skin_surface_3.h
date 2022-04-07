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

#ifndef CGAL_MESH_SKIN_SURFACE_3_H
#define CGAL_MESH_SKIN_SURFACE_3_H

#include <CGAL/license/Skin_surface_3.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/marching_tetrahedra_3.h>
#include <CGAL/Marching_tetrahedra_traits_skin_surface_3.h>
#include <CGAL/Marching_tetrahedra_observer_default_3.h>
#include <CGAL/Skin_surface_marching_tetrahedra_observer_3.h>
#include <CGAL/Skin_surface_polyhedral_items_3.h>

namespace CGAL {

template <class SkinSurface_3, class Polyhedron>
void mesh_skin_surface_3(SkinSurface_3 const &skin_surface, Polyhedron &p)
{
  skin_surface.mesh_surface_3(p);
}

//
//template <class SkinSurface_3, class P_Traits>
//void mesh_skin_surface_3
//(SkinSurface_3 const &skin_surface,
// Polyhedron_3<P_Traits, Skin_surface_polyhedral_items_3<SkinSurface_3> > &p)
//{
//  std::cout << "B" << std::endl;
//  skin_surface.mesh_skin_surface_3(p);
//}
//

} //namespace CGAL

#endif // CGAL_MESH_SKIN_SURFACE_3_H
