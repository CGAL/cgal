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

#ifndef CGAL_MESH_SKIN_SURFACE_3_H
#define CGAL_MESH_SKIN_SURFACE_3_H

#include <CGAL/Cartesian_converter.h>
#include <CGAL/marching_tetrahedra_3.h>
#include <CGAL/Marching_tetrahedra_traits_skin_surface_3.h>
#include <CGAL/Marching_tetrahedra_observer_default_3.h>
#include <CGAL/Skin_surface_marching_tetrahedra_observer_3.h>
#include <CGAL/Skin_surface_polyhedral_items_3.h>

CGAL_BEGIN_NAMESPACE

template <class SkinSurface_3, class Polyhedron>
void mesh_skin_surface_3(SkinSurface_3 const &skin_surface, 
			 Polyhedron &p)
{
  std::cout << "A" << std::endl;
  skin_surface.mesh_skin_surface_3(p);
//   typedef Marching_tetrahedra_traits_skin_surface_3<
//     SkinSurface_3,
//     typename SkinSurface_3::CMCT_Vertex_iterator,
//     typename SkinSurface_3::CMCT_Cell_iterator,
//     typename Polyhedron::HalfedgeDS>               Marching_tetrahedra_traits;
//   typedef Marching_tetrahedra_observer_default_3<
//     typename SkinSurface_3::CMCT_Vertex_iterator,
//     typename SkinSurface_3::CMCT_Cell_iterator,
//     Polyhedron>                                    Marching_tetrahedra_observer;

//   // Extract the coarse mesh using marching_tetrahedra
//   Marching_tetrahedra_traits   marching_traits(skin_surface);
//   Marching_tetrahedra_observer marching_observer;
//   marching_tetrahedra_3(skin_surface.cmct_vertices_begin(), 
// 			skin_surface.cmct_vertices_end(), 
// 			skin_surface.cmct_cells_begin(), 
// 			skin_surface.cmct_cells_end(), 
// 			p, 
// 			marching_traits,
// 			marching_observer);
}


template <class P_Traits,
	  class SkinSurface_3>
void mesh_skin_surface_3
(SkinSurface_3 const &skin_surface, 
 Polyhedron_3<P_Traits, 
 Skin_surface_polyhedral_items_3<SkinSurface_3> > &p)
{
  skin_surface.mesh_skin_surface_3(p);
//   typedef Polyhedron_3<P_Traits, 
//     Skin_surface_polyhedral_items_3<SkinSurface_3> > Polyhedron;

//   typedef Marching_tetrahedra_traits_skin_surface_3<
//     SkinSurface_3,
//     typename SkinSurface_3::CMCT_Vertex_iterator,
//     typename SkinSurface_3::CMCT_Cell_iterator,
//     typename Polyhedron::HalfedgeDS>               Marching_tetrahedra_traits;
//   typedef Marching_tetrahedra_observer_skin_surface_3<
//     typename SkinSurface_3::CMCT_Vertex_iterator,
//     typename SkinSurface_3::CMCT_Cell_iterator,
//     Polyhedron>                                    Marching_tetrahedra_observer;

//   // Extract the coarse mesh using marching_tetrahedra
//   Marching_tetrahedra_traits   marching_traits(skin_surface);
//   Marching_tetrahedra_observer marching_observer;
//   marching_tetrahedra_3(skin_surface.cmct_vertices_begin(), 
// 			skin_surface.cmct_vertices_end(), 
// 			skin_surface.cmct_cells_begin(), 
// 			skin_surface.cmct_cells_end(), 
// 			p, 
// 			marching_traits,
// 			marching_observer);
}



CGAL_END_NAMESPACE

#endif // CGAL_MESH_SKIN_SURFACE_3_H
