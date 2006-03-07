// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
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
// $Source: 
// $Revision: 1.1 $ $Date: 2005/12/12 16:20:58 $
// $Name:  $
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_MAKE_SURFACE_MESH_H
#define CGAL_MAKE_SURFACE_MESH_H

#include <CGAL/Surface_mesher/Surface_mesher.h>
#include <CGAL/Surface_mesh_traits_generator_3.h>
#include <CGAL/Surface_mesh_cell_base_3.h>
#include <CGAL/Surface_mesh_vertex_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Surface_mesh_default_criteria_3.h>


#include <CGAL/iterator.h> // CGAL::inserter()

namespace CGAL {

  struct Non_manifold_tag {};
  struct Manifold_tag {};
  struct Manifold_with_boundary_tag {};

template <typename C2T3,
	  typename Surface,
	  typename Criteria,
	  typename Tag>
void make_surface_mesh(C2T3& c2t3,
                       Surface surface,
                       Criteria criteria,
                       Tag tag = Non_manifold_tag() ) 
{
  typedef typename CGAL::Surface_mesh_traits_generator_3<Surface>::type Traits;

  make_surface_mesh(c2t3, surface, Traits(), criteria, tag);  
}

template <typename C2T3,
	  typename SurfaceMeshTraits_3,
	  typename Criteria>
void make_surface_mesh(C2T3& c2t3,
                       typename SurfaceMeshTraits_3::Surface_3 surface,
		       SurfaceMeshTraits_3 traits,
                       Criteria criteria,
                       Non_manifold_tag)
{
  using CGAL::Surface_mesher::Surface_mesher;
  
  typedef Surface_mesher<
    C2T3,
    typename SurfaceMeshTraits_3::Surface_3,
    SurfaceMeshTraits_3,
    Criteria> Mesher;

  typename SurfaceMeshTraits_3::Construct_initial_points get_initial_points =
    traits.construct_initial_points_object();

  get_initial_points(surface, CGAL::inserter(c2t3.triangulation()));
  Mesher mesher(c2t3, surface, traits, criteria);
  // TODO initial, then refine()
}

} // end namespace CGAL

#endif // CGAL_MAKE_SURFACE_MESH_H
